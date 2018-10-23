module PuMaSTutorials

using Weave

tutorial_directory = joinpath(@__DIR__,"..","tutorials")

function weave_all()
  tmp = joinpath(tutorial_directory,"multiple_response/multiple_response.jmd")
  #weave(tmp,doctype="pandoc",out_path=:pwd)
  weave(tmp,doctype = "md2html",out_path=:pwd)
  weave(tmp,doctype="md2pdf",out_path=:pwd)
end

end
