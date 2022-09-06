### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 50e81c16-2cd3-11ed-3f6e-c5ab57fc9565
begin
		using Pkg

	    Pkg.activate(pwd())
		Pkg.status()
end

# ╔═╡ 62f0ac59-db1c-4129-98ea-7b0fcc18a682
begin
	using PlutoUI
	PlutoUI.TableOfContents()
end

# ╔═╡ f0941b45-c60d-40fb-88cb-7bff02be8294
using MRIReco, BartIO, Plots

# ╔═╡ efc8c243-9dfa-4e2a-ad28-1414705cf06e
include("utils/utils_MP2RAGE.jl");

# ╔═╡ 4476fb51-2839-4e20-94db-5e5052802e39
md"# Prepare path"

# ╔═╡ a2bbaf22-07f5-40f0-ab54-34004125325a
pwd()

# ╔═╡ 53a7d708-e02d-44c9-bdf0-9b53850db592
if Sys.isapple() # my computer
	bart = BartIO.wrapper_bart("/Users/aurelien/Documents/Dev/mriSoft/bart")
else
	bart = BartIO.wrapper_bart("/home/runner/work/Benchmark_MRIReco_BART.jl/Benchmark_MRIReco_BART.jl/bart")
end

# ╔═╡ 721ba2e2-092b-4cd6-a6ca-50835ae97736
md"""# MP2RAGE
## Reconstruction of TI1/TI2"""

# ╔═╡ 7cb4abfc-647c-4454-8ab1-67edebc46e60
b = BrukerFile("data/MP2RAGE_FULLY/")

# ╔═╡ 86747c1b-7c23-45c3-a82f-07b461777f72
begin
	raw = acqDataFromMP2RAGE(b)
	acq = AcquisitionData(raw,OffsetBruker=true)
end

# ╔═╡ 42d02a7d-29d0-40ee-9640-a8db386144f4
begin
	p = Dict{Symbol, Any}()
	p[:reco] = "direct"
	p[:reconSize] = tuple(acq.encodingSize...)
	Ireco = reconstruction(acq, p)
	Isos = mergeChannels(Ireco)
end

# ╔═╡ e123ba17-b200-41e5-9eaf-521a878b947a
md"""## Estimate coil sensitivities
"""

# ╔═╡ b835fdc1-66ef-4f56-98b7-6788fa833d0f
begin
	kdata = extract3DKSpace(acq)
	kdata = permutedims(kdata,(1,2,3,5,4))
	sens = bart(1,"ecalib -m1 -c0 -r24",kdata[:,:,:,:,1])
end

# ╔═╡ 641c3e2a-49e7-469d-9f34-ce4f45ce9cf3
md"## Combine data with sensitivity and create MP2RAGE image"

# ╔═╡ 2ccfc121-7130-4e52-8bd9-8c0f70ee032b
size(Ireco)

# ╔═╡ 4950706d-313a-45a7-a9ae-ce1fac78ae80
begin
	sens_j = reshape(sens,acq.encodingSize[1],acq.encodingSize[2],acq.encodingSize[3],1,4);
	Isense = sum(conj.(sens_j).*Ireco,dims=5) ./sum(abs.(sens_j).^2,dims = 5)
	MP2_sense = mp2rage(Isense)
end

# ╔═╡ 87fd53b6-16de-42f6-9574-949f891a2a94
	p1 = heatmap( (MP2_sense[:,:,60,1,1]), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,title="fully",showaxis = false)

# ╔═╡ 8e11bcbf-4e55-47ca-959d-e4c9d64d5c10
md"# CS reconstruction
## Generate undersampling pattern"

# ╔═╡ bda923bf-9d72-4127-acd9-3190eedd5000
begin
	sx,sy,sz = acq.encodingSize
	maskBart = abs.(bart(1,"poisson -v -V5 -Y $sx -Z $sy -y1 -z1 -C 20"))
	maskBart = repeat(maskBart[1,:,:],1,1,sz)
	heatmap(maskBart[:,:,1])
end

# ╔═╡ 91153bab-8462-49be-a365-5c214415f810
begin
	# find indices
	I = findall(x->x==1,abs.(maskBart))
	subsampleInd = LinearIndices((sx,sy,sz))[I]
end;

# ╔═╡ 880b132b-7ae7-47c7-970e-ae5f1d9e9648
begin
	acqCS = deepcopy(acq);
	acqCS.subsampleIndices[1]=subsampleInd
	acqCS.subsampleIndices[2]=subsampleInd
	acqCS.kdata[1,1,1] = acqCS.kdata[1,1,1][subsampleInd,:]
	acqCS.kdata[2,1,1] = acqCS.kdata[2,1,1][subsampleInd,:]
end;

# ╔═╡ 56201036-cd4d-4b44-817f-3229405b8720
md"## Undersampling fft reco"

# ╔═╡ 6c1403a8-33d6-4424-a7c0-bdacb842ed06
begin
	Ireco_u = reconstruction(acqCS, p)
	Isos_u = mergeChannels(Ireco_u)
	Isense_u = sum(conj.(sens_j).*Ireco_u,dims=5) ./sum(abs.(sens_j).^2,dims = 5)
	MP2_u = mp2rage(Isense_u)
end

# ╔═╡ 9e502565-a5d7-489c-8cfc-f800be50c0b6
	p2 = heatmap( (MP2_u[:,:,60,1,1]), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,showaxis = false,title="undersampling")

# ╔═╡ 2ebf2b70-0566-47cd-a0ef-c591591b3963
md"## BART pics reconstruction"

# ╔═╡ 4240a7d3-559a-4f15-88ce-f87c24a975e8
begin
	kdata_u = extract3DKSpace(acqCS)
	kdata_u = permutedims(kdata_u,[1 2 3 5 4]);
	kdata_u = reshape(kdata_u,size(kdata_u)[1:end-1]...,1,size(kdata_u)[end])
end;

# ╔═╡ 7e98ff3a-08d7-421a-aac1-a47cfa1aaa49
heatmap(abs.(kdata_u[:,:,48,1,1,2]),clim = (0,5000))

# ╔═╡ 14200e12-21b6-405c-88a3-9c29f58da6c5
im_pics = bart(1,"pics -i 30 -S -R W:7:0:0.01",kdata_u,sens);

# ╔═╡ f697b417-a48e-495d-a6e4-4a7b29558f0a
MP2_pics = mp2rage(im_pics[:,:,:,1,1,:]);

# ╔═╡ 7ea69fa6-4c90-4b09-9114-4624771e32f4
begin
	p3 = heatmap( (MP2_pics[:,:,60,1,1]), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,showaxis = false, title="BART pics")
	plot(p1,p2,p3,layouts=(1,3))
end

# ╔═╡ 37b8bd03-c775-4eb7-b548-d8db43c04464
md"## Julia Fista reconstruction"

# ╔═╡ c6d3e734-ca53-41c5-bc56-0f0f98cdbb2e


# ╔═╡ Cell order:
# ╟─4476fb51-2839-4e20-94db-5e5052802e39
# ╠═62f0ac59-db1c-4129-98ea-7b0fcc18a682
# ╠═a2bbaf22-07f5-40f0-ab54-34004125325a
# ╠═50e81c16-2cd3-11ed-3f6e-c5ab57fc9565
# ╠═f0941b45-c60d-40fb-88cb-7bff02be8294
# ╠═efc8c243-9dfa-4e2a-ad28-1414705cf06e
# ╠═53a7d708-e02d-44c9-bdf0-9b53850db592
# ╟─721ba2e2-092b-4cd6-a6ca-50835ae97736
# ╠═7cb4abfc-647c-4454-8ab1-67edebc46e60
# ╠═86747c1b-7c23-45c3-a82f-07b461777f72
# ╠═42d02a7d-29d0-40ee-9640-a8db386144f4
# ╟─e123ba17-b200-41e5-9eaf-521a878b947a
# ╠═b835fdc1-66ef-4f56-98b7-6788fa833d0f
# ╟─641c3e2a-49e7-469d-9f34-ce4f45ce9cf3
# ╠═2ccfc121-7130-4e52-8bd9-8c0f70ee032b
# ╠═4950706d-313a-45a7-a9ae-ce1fac78ae80
# ╠═87fd53b6-16de-42f6-9574-949f891a2a94
# ╟─8e11bcbf-4e55-47ca-959d-e4c9d64d5c10
# ╠═bda923bf-9d72-4127-acd9-3190eedd5000
# ╠═91153bab-8462-49be-a365-5c214415f810
# ╠═880b132b-7ae7-47c7-970e-ae5f1d9e9648
# ╟─56201036-cd4d-4b44-817f-3229405b8720
# ╠═6c1403a8-33d6-4424-a7c0-bdacb842ed06
# ╠═9e502565-a5d7-489c-8cfc-f800be50c0b6
# ╟─2ebf2b70-0566-47cd-a0ef-c591591b3963
# ╠═4240a7d3-559a-4f15-88ce-f87c24a975e8
# ╠═7e98ff3a-08d7-421a-aac1-a47cfa1aaa49
# ╠═14200e12-21b6-405c-88a3-9c29f58da6c5
# ╠═f697b417-a48e-495d-a6e4-4a7b29558f0a
# ╠═7ea69fa6-4c90-4b09-9114-4624771e32f4
# ╟─37b8bd03-c775-4eb7-b548-d8db43c04464
# ╠═c6d3e734-ca53-41c5-bc56-0f0f98cdbb2e
