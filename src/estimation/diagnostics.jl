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
  [wresiduals(fpm, subjects[i], vvrandeffsorth[i], approx, fpm.args...; nsim=nsim, fpm.kwargs...) for i = 1:length(subjects)]
end
function wresiduals(fpm::FittedPumasModel, subject::Subject, randeffs, approx::LikelihoodApproximation, args...; nsim=nothing, kwargs...)
  is_sim = nsim == nothing
  if nsim == nothing
    approx = approx
    wres = wresiduals(fpm.model, subject, coef(fpm), randeffs, approx, args...; kwargs...)
    iwres = iwresiduals(fpm.model, subject, coef(fpm), randeffs, approx, args...; kwargs...)
  else
    approx = nothing
    wres = nothing
    iwres = eiwres(fpm.model, subject, coef(fpm), nsim)
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

"""
    restype(approx)

Returns the residual type for the given approximation method.
Can be one of [`FO`](@ref), [`FOCE`](@ref), or [`FOCEI`](@ref).
"""
restype(::FO) = :wres
restype(::FOCE) = :cwres
restype(::FOCEI) = :cwresi

function wresiduals(model::PumasModel, subject::Subject, param::NamedTuple, randeffs, approx::FO, args...; kwargs...)
  wres(model, subject, param, randeffs, args...; kwargs...)
end
function wresiduals(model::PumasModel, subject::Subject, param::NamedTuple, randeffs, approx::FOCE, args...; kwargs...)
  cwres(model, subject, param, randeffs, args...; kwargs...)
end
function wresiduals(model::PumasModel, subject::Subject, param::NamedTuple, randeffs, approx::FOCEI, args...; kwargs...)
  cwresi(model, subject, param, randeffs, args...; kwargs...)
end
function iwresiduals(model::PumasModel, subject::Subject, param::NamedTuple, randeffs, approx::FO, args...; kwargs...)
  iwres(model, subject, param, randeffs, args...; kwargs...)
end
function iwresiduals(model::PumasModel, subject::Subject, param::NamedTuple, randeffs, approx::FOCE, args...; kwargs...)
  icwres(model, subject, param, randeffs, args...; kwargs...)
end
function iwresiduals(model::PumasModel, subject::Subject, param::NamedTuple, randeffs, approx::FOCEI, args...; kwargs...)
  icwresi(model, subject, param, randeffs, args...; kwargs...)
end

"""

  wres(model, subject, param[, rfx])

To calculate the Weighted Residuals (WRES).
"""
function wres(m::PumasModel,
              subject::Subject,
              param::NamedTuple,
              vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
              args...;
              kwargs...)

  if vrandeffsorth isa Nothing
    vrandeffsorth=_orth_empirical_bayes(m, subject, param, FO(), args...; kwargs...)::AbstractVector
  end

  randeffstransform = totransform(m.random(param))
  randeffs = TransformVariables.transform(randeffstransform, vrandeffsorth)
  dist = _derived(m, subject, param, randeffs)

  F = _mean_derived_vηorth_jacobian(m, subject, param, vrandeffsorth, args...; kwargs...)

  _dv_keys = keys(subject.observations)
  return map(NamedTuple{_dv_keys}(_dv_keys)) do name

        V = Symmetric(F[name]*F[name]' + Diagonal(var.(dist[name])))
        return cholesky(V).U'\residuals(subject, dist)[name]
      end
end

"""
  cwres(model, subject, param[, rfx])

To calculate the Conditional Weighted Residuals (CWRES).
"""
function cwres(m::PumasModel,
               subject::Subject,
               param::NamedTuple,
               vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
               args...;
               kwargs...)

  if vrandeffsorth isa Nothing
    vrandeffsorth = _orth_empirical_bayes(m, subject, param, FOCE(), args...; kwargs...)
  end

  randeffstransform = totransform(m.random(param))
  randeffsEBE = TransformVariables.transform(randeffstransform, vrandeffsorth)

  distEBE = _derived(m, subject, param, randeffsEBE)

  _dv_keys = keys(subject.observations)
  foreach(_dv_keys) do _key
    if !_is_homoscedastic(distEBE[_key])
      throw(ArgumentError("dispersion parameter is not allowed to depend on the random effects when using FOCE"))
    end
    nothing
  end

  F = _mean_derived_vηorth_jacobian(m, subject, param, vrandeffsorth, args...; kwargs...)

  randeffstransform = totransform(m.random(param))

  return map(NamedTuple{_dv_keys}(_dv_keys)) do name
          V = Symmetric(F[name]*F[name]' + Diagonal(var.(distEBE[name])))
          return cholesky(V).U'\(residuals(subject, distEBE)[name] .+ F[name]*vrandeffsorth)
        end
end
"""
  cwresi(model, subject, param[, rfx])

To calculate the Conditional Weighted Residuals with Interaction (CWRESI).
"""

function cwresi(m::PumasModel,
                subject::Subject,
                param::NamedTuple,
                vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
                args...;
                kwargs...)

  if vrandeffsorth isa Nothing
    vrandeffsorth = _orth_empirical_bayes(m, subject, param, FOCEI(), args...; kwargs...)
  end

  randeffstransform = totransform(m.random(param))
  randeffs = TransformVariables.transform(randeffstransform, vrandeffsorth)
  dist = _derived(m, subject, param, randeffs)

  F = _mean_derived_vηorth_jacobian(m, subject, param, vrandeffsorth, args...; kwargs...)

  _dv_keys = keys(subject.observations)
  return map(NamedTuple{_dv_keys}(_dv_keys)) do name
           V = Symmetric(F[name]*F[name]' + Diagonal(var.(dist[name])))
           return cholesky(V).U'\(residuals(subject, dist)[name] .+ F[name]*vrandeffsorth)
         end
end

"""
  pred(model, subject, param[, rfx])

To calculate the Population Predictions (PRED).
"""
function pred(m::PumasModel,
              subject::Subject,
              param::NamedTuple,
              vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
              args...;
              kwargs...)

  if vrandeffsorth isa Nothing
    vrandeffsorth = _orth_empirical_bayes(m, subject, param, FO(), args...; kwargs...)
  end

  randeffs = TransformVariables.transform(totransform(m.random(param)), vrandeffsorth)
  dist = _derived(m, subject, param, randeffs)
  return map(d -> mean.(d), NamedTuple{keys(subject.observations)}(dist))
end


"""
  cpred(model, subject, param[, rfx])

To calculate the Conditional Population Predictions (CPRED).
"""
function cpred(m::PumasModel,
               subject::Subject,
               param::NamedTuple,
               vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
               args...;
               kwargs...)

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
          mean.(_derived(m, subject, param, _randeffs)[name])
        end,
        vrandeffsorth
      )
      return mean.(dist[name]) .- F*vrandeffsorth
    end
end

"""
  cpredi(model, subject, param[, rfx])

To calculate the Conditional Population Predictions with Interaction (CPREDI).
"""
function cpredi(m::PumasModel,
                subject::Subject,
                param::NamedTuple,
                vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
                args...;
                kwargs...)

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
  epred(model, subject, param, simulations_count)

To calculate the Expected Simulation based Population Predictions.
"""
function epred(m::PumasModel,
               subject::Subject,
               param::NamedTuple,
               nsim::Integer)
  sims = [simobs(m, subject, param).observed for i in 1:nsim]
  _dv_keys = keys(subject.observations)
  return map(name -> mean(getproperty.(sims, name)), NamedTuple{_dv_keys}(_dv_keys))
end

"""
  iwres(model, subject, param[, rfx])

To calculate the Individual Weighted Residuals (IWRES).
"""
function iwres(m::PumasModel,
               subject::Subject,
               param::NamedTuple,
               vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
               args...;
               kwargs...)

   if vrandeffsorth isa Nothing
     vrandeffsorth = _orth_empirical_bayes(m, subject, param, FO(), args...; kwargs...)
   end

  randeffs = TransformVariables.transform(totransform(m.random(param)), vrandeffsorth)
  dist = _derived(m, subject, param, randeffs)

  _dv_keys = keys(subject.observations)
  _res = residuals(subject, dist)
  return map(name -> _res[name] ./ std.(dist[name]), NamedTuple{_dv_keys}(_dv_keys))
end

"""
  icwres(model, subject, param[, rfx])

To calculate the Individual Conditional Weighted Residuals (ICWRES).
"""
function icwres(m::PumasModel,
                subject::Subject,
                param::NamedTuple,
                vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
                args...;
                kwargs...)

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

"""
  icwresi(model, subject, param[, rfx])

To calculate the Individual Conditional Weighted Residuals with Interaction (ICWRESI).
"""
function icwresi(m::PumasModel,
                 subject::Subject,
                 param::NamedTuple,
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
      sims_sum .+= (obsdv .- mean.(dv))./std.(dv)
    end
    return sims_sum ./ nsim
  end
end

function ipred(m::PumasModel,
                subject::Subject,
                param::NamedTuple,
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

function cipred(m::PumasModel,
                subject::Subject,
                param::NamedTuple,
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

function cipredi(m::PumasModel,
                  subject::Subject,
                  param::NamedTuple,
                  vrandeffsorth::Union{Nothing, AbstractVector}=nothing,
                  args...;
                  kwargs...)

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
  _icwresi = [icwresi(m, subject, param, vvrandeffsorth) for (subject, vvrandeffsorth) in zip(data, vvrandeffsorth)]
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
  _icwres = [icwres(m, subject, param, vvrandeffsorth) for (subject, vvrandeffsorth) in zip(data, vvrandeffsorth)]
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

function StatsBase.predict(model::PumasModel, subject::Subject, param, approx, vvrandeffsorth=_orth_empirical_bayes(model, subject, param, approx))
  pred = _predict(model, subject, param, approx, vvrandeffsorth)
  ipred = _ipredict(model, subject, param, approx, vvrandeffsorth)
  SubjectPrediction(pred, ipred, subject, approx)
end
StatsBase.predict(fpm::FittedPumasModel, approx::LikelihoodApproximation; kwargs...) = predict(fpm, fpm.data, approx; kwargs...)
function StatsBase.predict(fpm::FittedPumasModel, subjects::Population=fpm.data, approx=fpm.approx; nsim=nothing, timegrid=false,  useEBEs=true)
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
    return map(i -> predict(fpm.model, subjects[i], coef(fpm), approx, fpm.vvrandeffsorth[i]), 1:length(subjects))
  end
end

function StatsBase.predict(fpm::FittedPumasModel,
                           subject::Subject,
                           approx::LikelihoodApproximation=fpm.approx,
                           vrandeffsorth = _orth_empirical_bayes(fpm.model, subject, coef(fpm), approx);
                           timegrid=false)
  # We have not yet implemented custom time grids
  !(timegrid==false) && error("Using custom time grids is not yet implemented.")

  predict(fpm.model, subject, coef(fpm), approx, vrandeffsorth)
end

function _predict(model, subject, param, approx::FO, vvrandeffsorth)
  pred(model, subject, param, vvrandeffsorth)
end
function _ipredict(model, subject, param, approx::FO, vvrandeffsorth)
  ipred(model, subject, param, vvrandeffsorth)
end

function _predict(model, subject, param, approx::Union{FOCE, Laplace}, vvrandeffsorth)
  cpred(model, subject, param, vvrandeffsorth)
end
function _ipredict(model, subject, param, approx::Union{FOCE, Laplace}, vvrandeffsorth)
  cipred(model, subject, param, vvrandeffsorth)
end

function _predict(model, subject, param, approx::Union{FOCEI, LaplaceI}, vvrandeffsorth)
  cpredi(model, subject, param, vvrandeffsorth)
end
function _ipredict(model, subject, param, approx::Union{FOCEI, LaplaceI}, vvrandeffsorth)
  cipredi(model, subject, param, vvrandeffsorth)
end

function epredict(fpm, subject, nsim::Integer)
  epred(fpm.model, subjects, coef(fpm), nsim)
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

struct SubjectEBES{T1, T2, T3}
  ebes::T1
  subject::T2
  approx::T3
end
function empirical_bayes(fpm::FittedPumasModel, approx=fpm.approx)
  subjects = fpm.data

  trf = totransform(fpm.model.random(coef(fpm)))

  if approx == fpm.approx
    ebes = fpm.vvrandeffsorth
    return [SubjectEBES(TransformVariables.transform(trf, e), s, approx) for (e, s) in zip(ebes, subjects)]
  else
    # re-estimate under approx
    return [SubjectEBES(TransformVariables.transform(trf, _orth_empirical_bayes(fpm.model, subject, coef(fpm), approx), subject, approx)) for subject in subjects]
  end
end

function DataFrames.DataFrame(vebes::Vector{<:SubjectEBES}; include_covariates=true)
  subjects = [ebes.subject for ebes in vebes]
  df = select!(DataFrame(subjects; include_covariates=include_covariates, include_dvs=false), Not(:evid))
  for i = 1:length(first(vebes).ebes)
    df[!,Symbol("ebe_$i")] .= vcat((fill(ebes.ebes[i], length(ebes.subject.time)) for ebes in vebes)...)
  end
  df[!,:ebes_approx] .= vcat((fill(ebes.approx, length(ebes.subject.time)) for ebes in vebes)...)

  df
end

struct FittedPumasModelInspection{T1, T2, T3, T4}
  o::T1
  pred::T2
  wres::T3
  ebes::T4
end
StatsBase.predict(insp::FittedPumasModelInspection) = insp.pred
StatsBase.predict(insp::FittedPumasModelInspection, args...; kwargs...) = predict(insp.o, args...; kwargs...)
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
  FittedPumasModelInspection(fpm, pred, res, ebes)
end
function DataFrames.DataFrame(i::FittedPumasModelInspection; include_covariates=true)
  pred_df = DataFrame(i.pred; include_covariates=include_covariates)
  res_df = select!(select!(DataFrame(i.wres; include_covariates=false), Not(:id)), Not(:time))
  ebes_df = select!(select!(DataFrame(i.ebes; include_covariates=false), Not(:id)), Not(:time))

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
function findinfluential(m::PumasModel,
                         data::Population,
                         param::NamedTuple,

                         approx::LikelihoodApproximation,
                         args...;
                         k=5, kwargs...)
  d = [deviance(m, subject, param, approx, args...; kwargs...) for subject in data]
  p = partialsortperm(d, 1:k, rev=true)
  return [(data[pᵢ].id, d[pᵢ]) for pᵢ in p]
end
