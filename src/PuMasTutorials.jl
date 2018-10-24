module PuMaSTutorials

using Weave: weave

function weave_all(file)
    jmd = joinpath("tutorials", file)
    formats = map(dt -> joinpath(dt, replace(file, ".jmd" => string(".", dt))),
                  intersect(readdir(), ("html", "pdf", "md", "rst", "tex")))
    foreach(format -> weave(jmd, out_path = format), formats)
    return nothing
end

foreach(weave_all,
        (file for file in readdir("tutorials") if endswith(file, ".jmd")))

end
