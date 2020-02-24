function StatsBase.residuals(fpm::FittedPumasModel)
  # Return the residuals
  return [residuals(fpm.model, subject, coef(fpm), vrandeffsorth, fpm.args...; fpm.kwargs...) for (subject, vrandeffsorth) in zip(fpm.data, fpm.vvrandeffsorth)]
end
function StatsBase.residuals(model::PumasModel, subject::Subject, param::NamedTuple, vrandeffs::AbstractArray, args...; kwargs...)
  rtrf = totransform(model.random(param))
  randeffs = TransformVariables.transform(rtrf, vrandeffs)
  # Calculated the dependent variable distribution objects
  dist = _derived(model, subject, param, randeffs, args...; kwargs...)
  # Return the residuals
  return residuals(subject, dist)
end
function StatsBase.residuals(subject::Subject, dist)
  # Return the residuals
  _keys = keys(subject.observations)
  return map(x->x[1] .- mean.(x[2]), NamedTuple{_keys}(zip(subject.observations, dist)))
end
"""
  npde(model, subject, param, simulations_count)

To calculate the Normalised Prediction Distribution Errors (NPDE).
"""
function npde(m::PumasModel,
              subject::Subject,
              param::NamedTuple,
              nsim::Integer)

  _names = keys(subject.observations)
  sims = [simobs(m, subject, param).observed for i in 1:nsim]

  return map(NamedTuple{_names}(_names)) do name
           y = subject.observations[name]
           ysims = getproperty.(sims, name)
           mean_y = mean(ysims)
           cov_y = Symmetric(cov(ysims))
           Fcov_y = cholesky(cov_y)
           y_decorr = Fcov_y.U'\(y .- mean_y)

           φ = mean(ysims) do y_l
             y_decorr_l = Fcov_y.U'\(y_l .- mean_y)
             Int.(y_decorr_l .< y_decorr)
           end

           return quantile.(Normal(), φ)
         end
end

struct SubjectResidual{T1, T2, T3, T4}
  wres::T1
  iwres::T2
  subject::T3
  approx::T4
end
function wresiduals(fpm::FittedPumasModel, approx::LikelihoodApproximation=fpm.approx; nsim=nothing)
  subjects = fpm.data
  if approx == fpm.approx
    vvrandeffsorth = fpm.vvrandeffsorth
  else
    # re-estimate under approx
    vvrandeffsorth = [_orth_empirical_bayes(fpm.model, subject, coef(fpm), approx, fpm.args...; fpm.kwargs...) for subject in subjects]
  end
  [wresiduals(fpm, subjects[i], approx, vvrandeffsorth[i], fpm.args...; nsim=nsim, fpm.kwargs...) for i = 1:length(subjects)]
end
function wresiduals(
  fpm::FittedPumasModel,
  subject::Subject,
  approx::LikelihoodApproximation,
  vrandeffsorth::AbstractVector;
  nsim=nothing)

  is_sim = nsim == nothing
  if nsim == nothing
    approx = approx
    wres = wresiduals(fpm.model, subject, coef(fpm), approx, vrandeffsorth, fpm.args...; fpm.kwargs...)
    iwres = iwresiduals(fpm.model, subject, coef(fpm), approx, vrandeffsorth, fpm.args...; fpm.kwargs...)
  else
    approx = nothing
    wres = nothing
    iwres = eiwres(fpm.model, subject, coef(fpm), nsim, fpm.args...; fpm.kwargs...)
  end

  SubjectResidual(wres, iwres, subject, approx)
end

function DataFrames.DataFrame(vresid::Vector{<:SubjectResidual}; include_covariates=true)
  subjects = [resid.subject for resid in vresid]
  df = select!(DataFrame(subjects; include_covariates=include_covariates, include_dvs=false), Not(:evid))

  _keys = keys(first(subjects).observations)
  for name in (_keys)
    df[!,Symbol(string(name)*"_wres")] .= vcat((resid.wres[name] for resid in vresid)...)
    df[!,Symbol(string(name)*"_iwres")] .= vcat((resid.iwres[name] for resid in vresid)...)
    df[!,:wres_approx] .= vcat((fill(resid.approx, length(resid.subject.time)) for resid in vresid)...)
  end
  df
end

function wresiduals(
  model::PumasModel,
  subject::Subject,
  param::NamedTuple,
  approx::Union{FO,FOCE,FOCEI},
  vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
  args...; kwargs...)

  if vrandeffsorth isa Nothing
    vrandeffsorth = _orth_empirical_bayes(model, subject, param, approx, args...; kwargs...)
  end

  randeffstransform = totransform(model.random(param))
  randeffs = TransformVariables.transform(randeffstransform, vrandeffsorth)
  dist = _derived(model, subject, param, randeffs)

  F   = _mean_derived_vηorth_jacobian(model, subject, param, vrandeffsorth, args...; kwargs...)
  res = residuals(subject, dist)

  _dv_keys = keys(subject.observations)

  # For FOCE, we don't allow the dispersion parameter to depend on the random effects
  if approx isa FOCE
    foreach(_dv_keys) do _key
      if !_is_homoscedastic(dist[_key])
        throw(ArgumentError("dispersion parameter is not allowed to depend on the random effects when using FOCE"))
      end
      nothing
    end
  end

  return map(NamedTuple{_dv_keys}(_dv_keys)) do name
    # We have to handle missing values explicitly to avoid that the variance
    # components associated with missing values influence the weighting
    missingmask = ismissing.(subject.observations[name])
    Fname   = F[name]
    resname = res[name]
    Fname[missingmask, :] .= 0
    resname[missingmask]  .= 0

    V = Symmetric(Fname*Fname' + Diagonal(var.(dist[name])))

    # if "conditional" mothods, there is a first order term in the mean
    if approx isa Union{FOCE, FOCEI}
      resname .+= Fname*vrandeffsorth
    end

    ldiv!(cholesky(V).U', resname)

    # This setindex! operation fails when there are no missings even though
    # the mask is empty
    if any(missingmask)
      resname[missingmask] .= missing
    end

    return resname
  end
end

# PRED like _predict
function _predict(
  m::PumasModel,
  subject::Subject,
  param::NamedTuple,
  approx::FO,
  vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
  args...; kwargs...)

  if vrandeffsorth isa Nothing
    vrandeffsorth = _orth_empirical_bayes(m, subject, param, FO(), args...; kwargs...)
  end

  randeffs = TransformVariables.transform(totransform(m.random(param)), vrandeffsorth)
  dist = _derived(m, subject, param, randeffs)
  return map(d -> mean.(d), NamedTuple{keys(subject.observations)}(dist))
end


# CPRED like
function _predict(
  m::PumasModel,
  subject::Subject,
  param::NamedTuple,
  approx::FOCE,
  vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
  args...; kwargs...)

  if vrandeffsorth isa Nothing
    vrandeffsorth = _orth_empirical_bayes(m, subject, param, FOCE(), args...; kwargs...)
  end

  randeffstransform = totransform(m.random(param))
  randeffs = TransformVariables.transform(randeffstransform, vrandeffsorth)
  dist = _derived(m, subject, param, randeffs)


  _dv_keys = keys(subject.observations)
  return map(NamedTuple{_dv_keys}(_dv_keys)) do name
    F = ForwardDiff.jacobian(
      _vrandeffs -> begin
        _randeffs = TransformVariables.transform(randeffstransform, _vrandeffs)
        return mean.(_derived(m, subject, param, _randeffs)[name])
      end,
      vrandeffsorth
    )
    return mean.(dist[name]) .- F*vrandeffsorth
  end
end

# CPREDI like
function _predict(
  m::PumasModel,
  subject::Subject,
  param::NamedTuple,
  approx::Union{FOCEI, LaplaceI},
  vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
  args...; kwargs...)

  if vrandeffsorth isa Nothing
    vrandeffsorth = _orth_empirical_bayes(m, subject, param, FOCEI(), args...; kwargs...)
  end

  randeffstransform = totransform(m.random(param))
  randeffs = TransformVariables.transform(randeffstransform, vrandeffsorth)
  dist = _derived(m, subject, param, randeffs)
  _dv_keys = keys(subject.observations)
  return map(NamedTuple{_dv_keys}(_dv_keys)) do name
    F = ForwardDiff.jacobian(
      _vrandeffs -> begin
        _randeffs = TransformVariables.transform(randeffstransform, _vrandeffs)
        mean.(_derived(m, subject, param, _randeffs)[name])
      end,
      vrandeffsorth
    )
    return mean.(dist[name]) .- F*vrandeffsorth
  end
end

"""
  epredict(model, subject, param, simulations_count[, args...; kwarg...])

To calculate the Expected Simulation based Population Predictions.
"""
function epredict(
  m::PumasModel,
  subject::Subject,
  param::NamedTuple,
  nsim::Integer,
  args...; kwargs...)

  sims = [simobs(m, subject, param, args...; kwargs...).observed for i in 1:nsim]
  _dv_keys = keys(subject.observations)
  return map(name -> mean(getproperty.(sims, name)), NamedTuple{_dv_keys}(_dv_keys))
end

# IWRES like
function iwresiduals(
  m::PumasModel,
  subject::Subject,
  param::NamedTuple,
  approx::FO,
  vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
  args...; kwargs...)

   if vrandeffsorth isa Nothing
     vrandeffsorth = _orth_empirical_bayes(m, subject, param, FO(), args...; kwargs...)
   end

  randeffs = TransformVariables.transform(totransform(m.random(param)), vrandeffsorth)
  dist = _derived(m, subject, param, randeffs)

  _dv_keys = keys(subject.observations)
  _res = residuals(subject, dist)
  return map(name -> _res[name] ./ std.(dist[name]), NamedTuple{_dv_keys}(_dv_keys))
end

# ICWRES like
function iwresiduals(
  m::PumasModel,
  subject::Subject,
  param::NamedTuple,
  approx::FOCE,
  vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
  args...; kwargs...)

  if vrandeffsorth isa Nothing
    vrandeffsorth = _orth_empirical_bayes(m, subject, param, FOCE(), args...; kwargs...)
  end

  randeffstransform = totransform(m.random(param))
  randeffsEBE = TransformVariables.transform(randeffstransform, vrandeffsorth)
  dist = _derived(m, subject, param, randeffsEBE)

  _dv_keys = keys(subject.observations)

  foreach(_dv_keys) do _key
    if !_is_homoscedastic(dist[_key])
      throw(ArgumentError("dispersion parameter is not allowed to depend on the random effects when using FOCE"))
    end
    nothing
  end

  _res = residuals(subject, dist)
  return map(name -> _res[name] ./ std.(dist[name]), NamedTuple{_dv_keys}(_dv_keys))
end

# ICWRESI like
function iwresiduals(
  m::PumasModel,
  subject::Subject,
  param::NamedTuple,
  approx::Union{FOCEI,LaplaceI},
  vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
  args...;
  kwargs...)

  if vrandeffsorth isa Nothing
    vrandeffsorth = _orth_empirical_bayes(m, subject, param, FOCEI(), args...; kwargs...)
  end

  randeffs = TransformVariables.transform(totransform(m.random(param)), vrandeffsorth)
  dist = _derived(m, subject, param, randeffs)
  _dv_keys = keys(subject.observations)
  _res = residuals(subject, dist)
  return map(name -> _res[name] ./ std.(dist[name]), NamedTuple{_dv_keys}(_dv_keys))
end

"""
  eiwres(model, subject, param, simulations_count)

To calculate the Expected Simulation based Individual Weighted Residuals (EIWRES).
"""
function eiwres(m::PumasModel,
                subject::Subject,
                param::NamedTuple,
                nsim::Integer,
                args...;
                kwargs...)
  dist = _derived(m, subject, param, sample_randeffs(m, param), args...; kwargs...)
  _keys_dv = keys(subject.observations)
  return map(NamedTuple{_keys_dv}(_keys_dv)) do name
    dv = dist[name]
    obsdv = subject.observations[name]
    sims_sum = (obsdv .- mean.(dv))./std.(dv)
    for i in 2:nsim
      dist = _derived(m, subject, param, sample_randeffs(m, param), args...; kwargs...)
      dv = dist[name]
      sims_sum .+= (obsdv .- mean.(dv))./std.(dv)
    end
    return sims_sum ./ nsim
  end
end

function _ipredict(
  m::PumasModel,
  subject::Subject,
  param::NamedTuple,
  approx::FO,
  vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
  args...;
  kwargs...)

  if vrandeffsorth isa Nothing
    vrandeffsorth = _orth_empirical_bayes(m, subject, param, FO(), args...; kwargs...)
  end

  randeffs = TransformVariables.transform(totransform(m.random(param)), vrandeffsorth)
  dist = _derived(m, subject, param, randeffs, args...; kwargs...)
  return map(d->mean.(d), dist)
end

function _ipredict(
  m::PumasModel,
  subject::Subject,
  param::NamedTuple,
  approx::FOCE,
  vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
  args...;
  kwargs...)

  if vrandeffsorth isa Nothing
    vrandeffsorth = _orth_empirical_bayes(m, subject, param, FOCE(), args...; kwargs...)
  end

  randeffs = TransformVariables.transform(totransform(m.random(param)), vrandeffsorth)
  dist = _derived(m, subject, param, randeffs, args...; kwargs...)
  return map(d->mean.(d), dist)
end

function _ipredict(
  m::PumasModel,
  subject::Subject,
  param::NamedTuple,
  approx::Union{FOCEI, LaplaceI},
  vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
  args...; kwargs...)

  if vrandeffsorth isa Nothing
    vrandeffsorth = _orth_empirical_bayes(m, subject, param, FOCEI(), args...; kwargs...)
  end

  randeffs = TransformVariables.transform(totransform(m.random(param)), vrandeffsorth)
  dist = _derived(m, subject, param, randeffs, args...; kwargs...)
  return map(d->mean.(d), dist)
end

function ηshrinkage(m::PumasModel,
                    data::Population,
                    param::NamedTuple,
                    approx::LikelihoodApproximation,
                    args...;
                    kwargs...)

  vvrandeffsorth = [Pumas._orth_empirical_bayes(m, subject, param, approx, args...; kwargs...) for subject in data]
  vtrandeffs = [TransformVariables.transform(totransform(m.random(param)), _vrandefforth) for _vrandefforth in vvrandeffsorth]

  randeffsstd = map(keys(first(vtrandeffs))) do k
    return 1 .- std(getfield.(vtrandeffs, k)) ./ sqrt.(var(m.random(param).params[k]))
  end

  return NamedTuple{keys(first(vtrandeffs))}(randeffsstd)
end


function ϵshrinkage(m::PumasModel,
                    data::Population,
                    param::NamedTuple,
                    approx::FOCEI,
                    vvrandeffsorth::Union{Nothing, AbstractVector}=nothing,
                    args...;
                    kwargs...)

  if vvrandeffsorth isa Nothing
    vvrandeffsorth = [_orth_empirical_bayes(m, subject, param, FOCEI(), args...; kwargs...) for subject in data]
  end

  _keys_dv = keys(first(data).observations)
  _icwresi = [iwresiduals(m, subject, param, FOCEI(), vvrandeffsorth) for (subject, vvrandeffsorth) in zip(data, vvrandeffsorth)]
  map(name -> 1 - std(vec(VectorOfArray(getproperty.(_icwresi, name))), corrected = false), NamedTuple{_keys_dv}(_keys_dv))
end

function ϵshrinkage(m::PumasModel,
                    data::Population,
                    param::NamedTuple,
                    approx::FOCE,
                    vvrandeffsorth::Union{Nothing, AbstractVector}=nothing,
                    args...;
                    kwargs...)

  if vvrandeffsorth isa Nothing
    vvrandeffsorth = [_orth_empirical_bayes(m, subject, param, FOCE(), args...; kwargs...) for subject in data]
  end

  _keys_dv = keys(first(data).observations)
  _icwres = [iwresiduals(m, subject, param, FOCE(), vvrandeffsorth) for (subject, vvrandeffsorth) in zip(data, vvrandeffsorth)]
  map(name -> 1 - std(vec(VectorOfArray(getproperty.(_icwres, name))), corrected = false), NamedTuple{_keys_dv}(_keys_dv))
end

function StatsBase.aic(m::PumasModel,
                       data::Population,
                       param::NamedTuple,
                       approx::LikelihoodApproximation,
                       args...;
                       kwargs...)
  numparam = TransformVariables.dimension(totransform(m.param))
  2*(marginal_nll(m, data, param, approx, args...; kwargs...) + numparam)
end

function StatsBase.bic(m::PumasModel,
                       data::Population,
                       param::NamedTuple,
                       approx::LikelihoodApproximation,
                       args...;
                       kwargs...)
  numparam = TransformVariables.dimension(totransform(m.param))
  2*marginal_nll(m, data, param, approx, args...; kwargs...) + numparam*log(sum(t -> length(t.time), data))
end

### Predictions
struct SubjectPrediction{T1, T2, T3, T4}
  pred::T1
  ipred::T2
  subject::T3
  approx::T4
end

function StatsBase.predict(
  model::PumasModel,
  subject::Subject,
  param::NamedTuple,
  approx::LikelihoodApproximation,
  vvrandeffsorth::AbstractVector=_orth_empirical_bayes(model, subject, param, approx),
  args...; kwargs...
  )

  pred = _predict(model, subject, param, approx, vvrandeffsorth, args...; kwargs...)
  ipred = _ipredict(model, subject, param, approx, vvrandeffsorth, args...; kwargs...)
  SubjectPrediction(pred, ipred, subject, approx)
end

StatsBase.predict(fpm::FittedPumasModel, approx::LikelihoodApproximation; kwargs...) = predict(fpm, fpm.data, approx; kwargs...)

function StatsBase.predict(
  fpm::FittedPumasModel,
  subjects::Population=fpm.data,
  approx::LikelihoodApproximation=fpm.approx;
  nsim=nothing, timegrid=false, useEBEs=true)

  if !useEBEs
    error("Sampling from the omega distribution is not yet implemented.")
  end
  if !(timegrid==false)
    error("Using custom time grids is not yet implemented.")
  end
  if !(nsim isa Nothing)
    error("Using simulated subjects is not yet implemented.")
  end

  _estimate_bayes = approx == fpm.approx ? false : true

  if _estimate_bayes
    # re-estimate under approx
    return map(subject -> predict(fpm, subject, approx; timegrid=timegrid), subjects)
  else
    return map(i -> predict(fpm.model, subjects[i], coef(fpm), approx, fpm.vvrandeffsorth[i], fpm.args...; fpm.kwargs...), 1:length(subjects))
  end
end

function StatsBase.predict(
  fpm::FittedPumasModel,
  subject::Subject,
  approx::LikelihoodApproximation = fpm.approx,
  vrandeffsorth::AbstractVector = _orth_empirical_bayes(
    fpm.model,
    subject,
    coef(fpm),
    approx,
    fpm.args...; fpm.kwargs...);
  timegrid=false)
  # We have not yet implemented custom time grids
  !(timegrid==false) && error("Using custom time grids is not yet implemented.")

  predict(fpm.model, subject, coef(fpm), approx, vrandeffsorth, fpm.args...; fpm.kwargs...)
end

function DataFrames.DataFrame(vpred::Vector{<:SubjectPrediction}; include_covariates=true)
  subjects = [pred.subject for pred in vpred]
  df = select!(DataFrame(subjects; include_covariates=include_covariates, include_dvs=false), Not(:evid))
  _keys = keys(first(subjects).observations)
  for name in  _keys
    df[!,Symbol(string(name)*"_pred")] .= vcat((pred.pred[name] for pred in vpred)...)
    df[!,Symbol(string(name)*"_ipred")] .= vcat((pred.ipred[name] for pred in vpred)...)
    df[!,:pred_approx] .= vcat((fill(pred.approx, length(pred.subject.time)) for pred in vpred)...)
  end
  df
end

function empirical_bayes(fpm::FittedPumasModel, approx=fpm.approx)
  subjects = fpm.data

  trf = totransform(fpm.model.random(coef(fpm)))

  if approx == fpm.approx
    ebes = fpm.vvrandeffsorth
    return [TransformVariables.transform(trf, e) for (e, s) in zip(ebes, subjects)]
  else
    # re-estimate under approx
    return [
      TransformVariables.transform(
        trf,
        _orth_empirical_bayes(
          fpm.model, subject,
          coef(fpm),
          approx,
          fpm.args...; fpm.kwargs...
        )
      ) for subject in subjects]
  end
end

struct FittedPumasModelInspection{T1, T2, T3, T4, T5}
  o::T1
  pred::T2
  wres::T3
  ebes::T4
  ebes_approx::T5
end
StatsBase.predict(insp::FittedPumasModelInspection) = insp.pred
# We allow args... here since the called method will only use the saved args...
# from the fitting stage
StatsBase.predict(insp::FittedPumasModelInspection, args...) = predict(insp.o, args...)
wresiduals(insp::FittedPumasModelInspection) = insp.wres
empirical_bayes(insp::FittedPumasModelInspection) = insp.ebes

function inspect(fpm; pred_approx=fpm.approx, wres_approx=fpm.approx, ebes_approx=fpm.approx)
  print("Calculating: ")
  print("predictions")
  pred = predict(fpm, pred_approx)
  print(", weighted residuals")
  res = wresiduals(fpm, wres_approx)
  print(", empirical bayes")
  ebes = empirical_bayes(fpm, ebes_approx)
  println(". Done.")
  FittedPumasModelInspection(fpm, pred, res, ebes, ebes_approx)
end
function DataFrames.DataFrame(i::FittedPumasModelInspection; include_covariates=true)
  pred_df = DataFrame(i.pred; include_covariates=include_covariates)
  res_df = select!(select!(DataFrame(i.wres; include_covariates=false), Not(:id)), Not(:time))

  ebe_keys = keys(first(i.ebes))
  ebe_types = map(typeof, first(i.ebes))
  ebes_df = select!(DataFrame(i.o.data; include_covariates=false, include_dvs=false), Not(:evid))
  for k ∈ ebe_keys
    ebe_type = ebe_types[k]
    if ebe_type <: Number
      ebes_df[!, k] .= vcat((fill(i.ebes[n][k], length(i.o.data[n].time)) for n = 1:length(i.o.data))...)
    elseif ebe_type <: AbstractVector
      for j = 1:length(first(i.ebes)[k])
        ebes_df[!, Symbol(string(k)*"_$j")] .= vcat((fill(i.ebes[n][k][j], length(i.o.data[n].time)) for n = 1:length(i.o.data))...)
      end
    end
  end
  ebes_df[!, :ebes_approx] .= summary(i.ebes_approx)
  ebes_df = select!(select!(ebes_df, Not(:id)), Not(:time))
  df = hcat(pred_df, res_df, ebes_df)
end




################################################################################
#                              Plotting functions                              #
################################################################################

########################################
#   Convergence plot infrastructure    #
########################################

"""
    _objectivefunctionvalues(obj)

Returns the objective function values during optimization.
Must return a `Vector{Number}`.
"""
_objectivefunctionvalues(f::FittedPumasModel) = getproperty.(f.optim.trace, :value)

"""
    _convergencedata(obj; metakey="x")

Returns the "timeseries" of optimization as a matrix, with series as columns.
!!! warn
    This must return parameter data in the same order that [`_paramnames`](@ref)
    returns names.
"""
function _convergencedata(f::FittedPumasModel; metakey="x")

  metakey != "x" && return transpose(hcat(getindex.(getproperty.(f.optim.trace, :metadata), metakey)...))

  trf  = totransform(f.model.param)         # get the transform which has been applied to the params
  itrf = toidentitytransform(f.model.param) # invert the param transform

  return transpose(                                     # return series as columns
              hcat(TransformVariables.inverse.(         # apply the inverse of the given transform to the data.
                  Ref(itrf),                            # wrap in a `Ref`, to avoid broadcasting issues
                  TransformVariables.transform.(        # apply the initial transform to the process
                      Ref(trf),                         # again - make sure no broadcasting across the `TransformTuple`
                      getindex.(                        # get every `x` vector from the metadata of the trace
                          getproperty.(                 # get the metadata of each trace element
                              f.optim.trace, :metadata  # getproperty expects a `Symbol`
                              ),
                          metakey                           # property x is a key for a `Dict` - hence getindex
                          )
                      )
                  )...                                  # splat to get a matrix out
              )
          )
end

"""
    _paramnames(obj)

Returns the names of the parameters which convergence is being checked for.
!!! warn
    This must return parameter names in the same order that [`_convergencedata`](@ref)
    returns data.
"""
function _paramnames(f::FittedPumasModel)
  paramnames = [] # empty array, will fill later
  for (paramname, paramval) in pairs(coef(f)) # iterate through the parameters
    # decompose all parameters (matrices, etc.) into scalars and name them appropriately
    _push_varinfo!(paramnames, [], nothing, nothing, paramname, paramval, nothing, nothing)
  end
  return paramnames
end

# This will use default args and kwargs!!
findinfluential(fpm::FittedPumasModel) = findinfluential(fpm.model, fpm.data, coef(fpm), fpm.approx, fpm.args...; fpm.kwargs...)
function findinfluential(
  m::PumasModel,
  data::Population,
  param::NamedTuple,
  approx::LikelihoodApproximation,
  args...;
  k=5, kwargs...)

  d = [deviance(m, subject, param, approx, args...; kwargs...) for subject in data]
  p = partialsortperm(d, 1:k, rev=true)
  return [(data[pᵢ].id, d[pᵢ]) for pᵢ in p]
end
