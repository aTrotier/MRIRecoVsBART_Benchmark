### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 5c30fca3-c0af-4919-9978-a19be968f781
begin
		using Pkg
	
	    Pkg.activate(pwd())
		Pkg.instantiate()
		Pkg.status()	
end

# ╔═╡ a76cbc7a-0d73-4041-a7ea-18ca1e75cd84
begin
	using MRIReco, BartIO, Plots,PlutoUI
	PlutoUI.TableOfContents()
end

# ╔═╡ f2e90faa-2e1a-11ed-0709-d54f5d6ad389
pwd()

# ╔═╡ f0b62eb2-49a7-4b17-8af7-85986b750b10
if Sys.isapple()
	bart = BartIO.wrapper_bart("/Users/aurelien/Documents/Dev/mriSoft/bart")
else
	bart = BartIO.wrapper_bart("/home/runner/work/Benchmark_MRIReco_BART.jl/Benchmark_MRIReco_BART.jl/bart")
end

# ╔═╡ 682f31df-bf2e-4b1e-96f0-a4ab4606114a
md"# Test Bart in pluto notebook"

# ╔═╡ 612fe7fa-6fa8-48dd-a947-5866cc0ffa27
phantom = bart(1,"phantom")

# ╔═╡ ec957a2d-9810-46d9-87fb-bad4252b365f
heatmap(abs.(phantom), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,showaxis = false)

# ╔═╡ Cell order:
# ╠═f2e90faa-2e1a-11ed-0709-d54f5d6ad389
# ╠═5c30fca3-c0af-4919-9978-a19be968f781
# ╠═a76cbc7a-0d73-4041-a7ea-18ca1e75cd84
# ╠═f0b62eb2-49a7-4b17-8af7-85986b750b10
# ╟─682f31df-bf2e-4b1e-96f0-a4ab4606114a
# ╠═612fe7fa-6fa8-48dd-a947-5866cc0ffa27
# ╠═ec957a2d-9810-46d9-87fb-bad4252b365f
