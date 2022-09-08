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
	using ImageQualityIndexes:assess_ssim
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

# ╔═╡ 8e92b598-0ec7-4e45-879e-98badc57fee6
md"# BART processing"

# ╔═╡ fac005c7-4965-4e0a-a828-cfdd67aae4cf
md"## Generate 3D phantom"

# ╔═╡ 46534d70-14c2-4565-8d98-7ac937270c0f
begin
	sx = 128
	nch = 8
end

# ╔═╡ 1526b4f0-5a9a-4d82-9392-2fb97e62530c
begin
	phant3D = bart(1,"phantom -3 -x$sx -s$nch");
	phant3D = bart(1,"noise -n0.005",phant3D)
	phant3D_rss = bart(1,"rss 8",phant3D)
end;

# ╔═╡ b6cb859d-d9c6-4f28-ae9e-a04bc912b2d2
heatmap(abs.(phant3D_rss[:,:,80]), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,showaxis = false)

# ╔═╡ e62e9218-918d-4341-b26f-e4bd2df9d18d
kbart = bart(1,"fft -u 7",phant3D);

# ╔═╡ d73ab523-3840-4b33-a118-1bb0d650e306
md"## create mask and undersample"

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

# ╔═╡ 614efa5f-4ca0-441f-b519-527fb5564d92
md"## Estimates sensitivity maps"

# ╔═╡ 56bd19d5-f7e6-465b-9bf5-88a834efa2c4
sens = bart(1,"ecalib -m1 -c0",kbart_u);

# ╔═╡ 7428391b-95a1-47a6-93c6-dc8208fbac28
md"## PICS reconstruction"

# ╔═╡ 912de2f2-cdb3-43ef-9077-dff56d2f82d7
t1 = @elapsed im_pics = bart(1,"pics -d5 -i30 -RW:7:0:0.01",kbart_u,sens)

# ╔═╡ e4e8a554-2e3d-4479-901d-e5f7179eada4
	heatmap(abs.(im_pics[:,:,80]), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,showaxis = false)

# ╔═╡ 4db123bf-f4fb-492a-b987-40b63a95af70
imFully = bart(1,"pics -d5 -i1 l2 -r0",kbart,sens);

# ╔═╡ c255089a-23d7-4c42-895d-57f0cd840ac9
RMSE_bart = MRIReco.norm(vec(abs.(im_pics))-vec(abs.(imFully)))/MRIReco.norm(vec(abs.(imFully)))

# ╔═╡ bb1f0a12-c707-4733-8112-7eac58d14ddb
ssim_bart = round(assess_ssim(abs.(im_pics[:,:,80]),abs.(imFully[:,:,80])),digits=3)

# ╔═╡ 01e9bbef-ba51-4a09-9598-4ea7a1403764
md"# Julia processing"

# ╔═╡ b580f18f-03be-4e07-a916-5d8851466220
md"## Fully sampled datasets"

# ╔═╡ 5aaad47b-ef42-4dbe-9011-8829b4e89f7f
tr = MRIBase.CartesianTrajectory3D(Float64,128,128,numSlices=128)

# ╔═╡ 2d8a3d57-db3b-4cda-9a2b-9031ba6c5ede
kdata_j = [reshape(ComplexF64.(kbart),:,nch) for i=1:1, j=1:1, k=1:1]

# ╔═╡ 0bb2bc6d-8afc-41a9-87ec-c7f820cc02bd
begin
	acq = AcquisitionData(tr,kdata_j)
	acq.encodingSize = [sx,sx,sx]
end

# ╔═╡ c981d2c8-c81e-4072-8f96-91cc11615bef
begin
	params = Dict{Symbol, Any}()
	params[:reco] = "direct"
	params[:reconSize] = tuple(acq.encodingSize...)
	Ireco = reconstruction(acq, params)
	Isos = mergeChannels(Ireco)
	heatmap(abs.(Isos[:,:,80]), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,showaxis = false)
end

# ╔═╡ 32abd59d-5d66-4a0d-9f59-5207a3b88b75
md"## Undersampled the data"

# ╔═╡ 1aed83a6-33f8-451b-88cf-62744d905cbe
begin
	# find indices
	I = findall(x->x==1,abs.(repeat(mask,sx,1,1)))
	subsampleInd = LinearIndices((sx,sx,sx))[I]
end;

# ╔═╡ 9256ca72-0342-4e09-88d1-825622567501
begin
	acqCS = deepcopy(acq);
	acqCS.subsampleIndices[1]=subsampleInd
	acqCS.kdata[1,1,1] = acqCS.kdata[1,1,1][subsampleInd,:]
end;

# ╔═╡ aa168a9c-bc03-4b0c-ad3c-264fd4ebfc5a
begin
	Ireco_u = reconstruction(acqCS, params)
	Isos_u = mergeChannels(Ireco_u)
	heatmap(abs.(Isos_u[:,:,80]), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,showaxis = false)
end

# ╔═╡ 544f06aa-5029-45e0-858d-4cdf55754dd1
md"## Reco fista"

# ╔═╡ fa3fb104-5a4a-4841-b4ee-ba3ee45540b0
begin
	params2 = Dict{Symbol, Any}()
	params2[:reco] = "multiCoil"
	params2[:reconSize] = tuple(acqCS.encodingSize...)
	params2[:senseMaps] = ComplexF64.(sens);
	
	params2[:solver] = "fista"
	params2[:sparseTrafoName] = "Wavelet"
	params2[:regularization] = "L1"
	params2[:λ] = 0.01 # 5.e-2
	params2[:iterations] = 30
	params2[:normalize_ρ] = false
	params2[:ρ] = 0.95
	#params2[:relTol] = 0.1
	params2[:normalizeReg] = true
	
	
	t2 = @elapsed I_wav = reconstruction(acqCS, params2);
end

# ╔═╡ f02c089e-0d92-4984-a315-c378adea245f
	heatmap(abs.(I_wav[:,:,80]), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,showaxis = false)

# ╔═╡ 9b3f09d1-d1ae-4013-a057-a48a4cacaead
RMSE_julia = MRIReco.norm(vec(abs.(I_wav))-vec(abs.(Isos)))/MRIReco.norm(vec(abs.(Isos)))

# ╔═╡ 55579a0b-12e8-471a-891c-a1151df7342a
ssim_julia = round(assess_ssim(abs.(I_wav[:,:,80]),abs.(Isos[:,:,80])),digits=3)

# ╔═╡ 5cc0323c-c255-4176-87c9-2bdf306592af
md" # Benchmark Results"

# ╔═╡ abd6ac0d-123e-4463-97b5-877552a240a1
md"## Reconstruction time"

# ╔═╡ d57486a6-0175-4135-9737-cb610ad82a83
begin #round for plots
	t_bart = round(t1)
	t_julia = round(t2)
	rapport = round(t2/t1,digits=2)

	RMSE_bart2=round(Float64.(RMSE_bart),digits =4)
	RMSE_julia2=round(RMSE_julia,digits =4)
end;

# ╔═╡ 1f9a0b20-63b0-4198-a93e-114ba881571f
md"Reconstruction time julia = $t_julia sec vs t_bart = $t1 sec

**Bart is $rapport faster**"

# ╔═╡ da073e09-8375-496c-8d18-889d2d563eff
md"## Precision :

**RMSE:**

julia = $RMSE_julia2

bart = $RMSE_bart2

**SSIM:**

julia = $ssim_julia

bart = $ssim_bart"

# ╔═╡ Cell order:
# ╠═f2e90faa-2e1a-11ed-0709-d54f5d6ad389
# ╠═5c30fca3-c0af-4919-9978-a19be968f781
# ╠═a76cbc7a-0d73-4041-a7ea-18ca1e75cd84
# ╠═f0b62eb2-49a7-4b17-8af7-85986b750b10
# ╟─682f31df-bf2e-4b1e-96f0-a4ab4606114a
# ╠═612fe7fa-6fa8-48dd-a947-5866cc0ffa27
# ╠═ec957a2d-9810-46d9-87fb-bad4252b365f
# ╟─8e92b598-0ec7-4e45-879e-98badc57fee6
# ╟─fac005c7-4965-4e0a-a828-cfdd67aae4cf
# ╠═46534d70-14c2-4565-8d98-7ac937270c0f
# ╠═1526b4f0-5a9a-4d82-9392-2fb97e62530c
# ╠═b6cb859d-d9c6-4f28-ae9e-a04bc912b2d2
# ╠═e62e9218-918d-4341-b26f-e4bd2df9d18d
# ╟─d73ab523-3840-4b33-a118-1bb0d650e306
# ╠═c4490b95-527e-49a2-9b90-c31c86c82dc0
# ╠═c922fe4e-4a61-4a82-a9ae-1421b7c40516
# ╠═1493aea6-db34-44fb-b128-81ea92a061c8
# ╠═5daabe9e-6d1e-4b53-8f37-c42a56d431af
# ╟─614efa5f-4ca0-441f-b519-527fb5564d92
# ╠═56bd19d5-f7e6-465b-9bf5-88a834efa2c4
# ╟─7428391b-95a1-47a6-93c6-dc8208fbac28
# ╠═912de2f2-cdb3-43ef-9077-dff56d2f82d7
# ╠═e4e8a554-2e3d-4479-901d-e5f7179eada4
# ╠═4db123bf-f4fb-492a-b987-40b63a95af70
# ╠═c255089a-23d7-4c42-895d-57f0cd840ac9
# ╠═bb1f0a12-c707-4733-8112-7eac58d14ddb
# ╟─01e9bbef-ba51-4a09-9598-4ea7a1403764
# ╟─b580f18f-03be-4e07-a916-5d8851466220
# ╠═5aaad47b-ef42-4dbe-9011-8829b4e89f7f
# ╠═2d8a3d57-db3b-4cda-9a2b-9031ba6c5ede
# ╠═0bb2bc6d-8afc-41a9-87ec-c7f820cc02bd
# ╠═c981d2c8-c81e-4072-8f96-91cc11615bef
# ╟─32abd59d-5d66-4a0d-9f59-5207a3b88b75
# ╠═1aed83a6-33f8-451b-88cf-62744d905cbe
# ╠═9256ca72-0342-4e09-88d1-825622567501
# ╠═aa168a9c-bc03-4b0c-ad3c-264fd4ebfc5a
# ╠═544f06aa-5029-45e0-858d-4cdf55754dd1
# ╠═fa3fb104-5a4a-4841-b4ee-ba3ee45540b0
# ╠═f02c089e-0d92-4984-a315-c378adea245f
# ╠═9b3f09d1-d1ae-4013-a057-a48a4cacaead
# ╠═55579a0b-12e8-471a-891c-a1151df7342a
# ╟─5cc0323c-c255-4176-87c9-2bdf306592af
# ╟─abd6ac0d-123e-4463-97b5-877552a240a1
# ╠═d57486a6-0175-4135-9737-cb610ad82a83
# ╟─1f9a0b20-63b0-4198-a93e-114ba881571f
# ╟─da073e09-8375-496c-8d18-889d2d563eff
