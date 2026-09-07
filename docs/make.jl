using Documenter
using EnsembleMCMC

DocMeta.setdocmeta!(EnsembleMCMC, :DocTestSetup, :(using EnsembleMCMC); recursive=true)

makedocs(
    root = @__DIR__,
    sitename = "EnsembleMCMC.jl",
    modules = [EnsembleMCMC],
    checkdocs = :exports,
    doctest = true,
    warnonly = false,
    format = Documenter.HTML(prettyurls=true, edit_link="main"),
    pages = ["Getting started" => "index.md", "API" => "api.md"],
)
