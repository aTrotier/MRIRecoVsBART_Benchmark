### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 5c30fca3-c0af-4919-9978-a19be968f781
begin
		using Pkg
	
	    Pkg.activate(pwd())
		Pkg.instantiate()
end;

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
	bart = BartIO.wrapper_bart("/home/runner/work/MRIRecoVsBART_Benchmark/MRIRecoVsBART_Benchmark/bart")
end

# ╔═╡ 682f31df-bf2e-4b1e-96f0-a4ab4606114a
md"# Test Bart in pluto notebook"

# ╔═╡ 612fe7fa-6fa8-48dd-a947-5866cc0ffa27
begin
	phantom = bart(1,"phantom");
	phantom = bart(1,"noise -n0.005",phantom)
end

# ╔═╡ ec957a2d-9810-46d9-87fb-bad4252b365f
heatmap(abs.(phantom), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,showaxis = false)

# ╔═╡ fac005c7-4965-4e0a-a828-cfdd67aae4cf
md"# Generate 3D phantom"

# ╔═╡ 1526b4f0-5a9a-4d82-9392-2fb97e62530c
begin
	phant3D = bart(1,"phantom -3 -x128 -s8");
	phant3D = bart(1,"noise -n0.005",phant3D)
end;

# ╔═╡ b6cb859d-d9c6-4f28-ae9e-a04bc912b2d2
heatmap(abs.(phant3D[:,:,80,8]), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,showaxis = false)

# ╔═╡ e62e9218-918d-4341-b26f-e4bd2df9d18d
kbart = bart(1,"fft 7",phant3D);

# ╔═╡ c4490b95-527e-49a2-9b90-c31c86c82dc0
mask = bart(1,"poisson -Y 128 -Z 128 -y1.2 -z1.2 -C 20 -v -V5");

# ╔═╡ c922fe4e-4a61-4a82-a9ae-1421b7c40516
heatmap(abs.(mask[1,:,:]), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,showaxis = false)

# ╔═╡ 1493aea6-db34-44fb-b128-81ea92a061c8
kbart_u = kbart .* mask;

# ╔═╡ 5daabe9e-6d1e-4b53-8f37-c42a56d431af
begin
	im_u = bart(1,"fft -i 7",kbart_u)
	im_u_sos = bart(1,"rss 8",im_u)
	heatmap(abs.(im_u_sos[:,:,80]), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,showaxis = false)
end

# ╔═╡ 56bd19d5-f7e6-465b-9bf5-88a834efa2c4
sens = bart(1,"ecalib -m1 -c0",kbart_u);

# ╔═╡ 912de2f2-cdb3-43ef-9077-dff56d2f82d7
im_pics = bart(1,"pics -i30 -S -RW:7:0:0.01",kbart_u,sens);

# ╔═╡ e4e8a554-2e3d-4479-901d-e5f7179eada4
	heatmap(abs.(im_pics[:,:,80]), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,showaxis = false)

# ╔═╡ Cell order:
# ╠═f2e90faa-2e1a-11ed-0709-d54f5d6ad389
# ╠═5c30fca3-c0af-4919-9978-a19be968f781
# ╠═a76cbc7a-0d73-4041-a7ea-18ca1e75cd84
# ╠═f0b62eb2-49a7-4b17-8af7-85986b750b10
# ╟─682f31df-bf2e-4b1e-96f0-a4ab4606114a
# ╠═612fe7fa-6fa8-48dd-a947-5866cc0ffa27
# ╠═ec957a2d-9810-46d9-87fb-bad4252b365f
# ╟─fac005c7-4965-4e0a-a828-cfdd67aae4cf
# ╠═1526b4f0-5a9a-4d82-9392-2fb97e62530c
# ╠═b6cb859d-d9c6-4f28-ae9e-a04bc912b2d2
# ╠═e62e9218-918d-4341-b26f-e4bd2df9d18d
# ╠═c4490b95-527e-49a2-9b90-c31c86c82dc0
# ╠═c922fe4e-4a61-4a82-a9ae-1421b7c40516
# ╠═1493aea6-db34-44fb-b128-81ea92a061c8
# ╠═5daabe9e-6d1e-4b53-8f37-c42a56d431af
# ╠═56bd19d5-f7e6-465b-9bf5-88a834efa2c4
# ╠═912de2f2-cdb3-43ef-9077-dff56d2f82d7
# ╠═e4e8a554-2e3d-4479-901d-e5f7179eada4
