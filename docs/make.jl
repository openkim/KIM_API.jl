using Documenter, KIM_API

# Set up documentation
makedocs(
    modules = [KIM_API],
    sitename = "KIM_API.jl",
    authors = "Amit Gupta <gupta839@umn.edu>",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://openkim.github.io/KIM_API.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Molly.jl Integration" => "molly.md",
        "Adding KIM Support to Your Simulator" => "simulator_integration.md",
        "API Reference" => [
            "High-level Interface" => "api/highlevel.md",
            "Library Management" => "api/libkim.md",
            "Model Management" => "api/model.md",
            "Species Handling" => "api/species.md",
            "Neighbor Lists" => "api/neighborlist.md",
            "Constants & Units" => "api/constants.md",
            "Utilities" => "api/utils.md",
        ],
        "Examples" => "examples.md",
        # "Troubleshooting" => "troubleshooting.md",
        # "Developer Guide" => "developer.md"
    ],
    repo = "https://github.com/openkim/KIM_API.jl/blob/{commit}{path}#L{line}",
    clean = true,
    doctest = false,
    linkcheck = false,
)

# Deploy documentation
deploydocs(
    repo = "github.com/openkim/KIM_API.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "master",
)
