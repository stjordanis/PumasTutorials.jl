module PuMaSTutorials

using Weave: weave

function weave_all(file)
    jmd = joinpath("tutorials", file)
    weave(jmd, doctype = "md2html", out_path = "html")
    weave(jmd, doctype = "md2pdf", out_path = "pdf")
    return nothing
end

foreach(weave_all,
        file for file in readdir("tutorials") if endswith(file, ".jmd"))

end
