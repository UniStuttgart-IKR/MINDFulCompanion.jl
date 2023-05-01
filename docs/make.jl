using Documenter, MINDFulCompanion

makedocs(sitename="MINDFulCompanion.jl",
    pages = [
        "Introduction" => "index.md",
        "API" => "API.md"
    ])

 deploydocs(
     repo = "github.com/UniStuttgart-IKR/MINDFulCompanion.jl.git",
 )
