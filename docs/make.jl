using Documenter, FinEtools, FinEtoolsMultithreading

makedocs(
	modules = [FinEtoolsMultithreading],
	doctest = false, clean = true,
	warnonly = Documenter.except(:linkcheck, :footnote),
	format = Documenter.HTML(prettyurls = false),
	authors = "Petr Krysl",
	sitename = "FinEtoolsMultithreading.jl",
	pages = Any[
	"Home" => "index.md",
	"How to guide" => "guide/guide.md",
	"Types and Functions" => Any[
		"man/man.md"]
		]
	)

deploydocs(
    repo = "github.com/PetrKryslUCSD/FinEtoolsMultithreading.jl.git",
)
