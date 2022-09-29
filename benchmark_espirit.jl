### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ e979c3be-3f37-11ed-13f0-d966a2fe1aba
begin
	using Pkg

	    Pkg.activate(pwd())
		Pkg.instantiate()
end;

# ╔═╡ b9cf2100-3959-454a-9ba5-d49bbb0c826c
begin
	using MRIReco, BartIO, Plots,PlutoUI
	PlutoUI.TableOfContents()
end

# ╔═╡ b25c25f1-3154-40f5-a58d-f945a29ed5f8
begin
    display(Threads.nthreads())
    ENV["OMP_NUM_THREADS"]=2 #Threads.nthreads() when FFTop will use them
    if Sys.isapple()
        bart = BartIO.wrapper_bart("/Users/aurelien/Documents/Dev/mriSoft/bart")
    else
        bart = BartIO.wrapper_bart("/home/runner/work/MRIRecoVsBART_Benchmark/MRIRecoVsBART_Benchmark/bart")
    end
end

# ╔═╡ 8087d5a8-99eb-4c11-97f2-883e41bce5b6
begin
	sx = 128
	nch = 8
	
	phant3D = bart(1,"phantom -3 -x$sx -s$nch");
	phant3D = bart(1,"noise -n0.005",phant3D)
	phant3D_rss = bart(1,"rss 8",phant3D)
	kbart = bart(1,"fft -u 7",phant3D);
end;

# ╔═╡ 550c0b48-0d94-4bf2-a15b-0cebcffe58bb
heatmap(abs.(phant3D_rss[:,:,80]), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,showaxis = false)

# ╔═╡ 3951341a-f2e7-465e-b783-a92d56bc05e1
t_bart = @elapsed sens_bart = bart(1,"ecalib -m1 -c0 -k 6 -r24",kbart);

# ╔═╡ c7295363-31e2-4053-8d8c-abeed4274a89
"""
Documentation : Crop the central area for 4D array
"""
function crop(A::Array{T,4}, s::NTuple{3,Int64}) where {T}
    nx, ny, nz = size(A)
    idx_x = div(nx, 2)-div(s[1], 2)+1:div(nx, 2)-div(s[1], 2)+s[1]
    idx_y = div(ny, 2)-div(s[2], 2)+1:div(ny, 2)-div(s[2], 2)+s[2]
    idx_z = div(nz, 2)-div(s[3], 2)+1:div(nz, 2)-div(s[3], 2)+s[3]
    return A[idx_x, idx_y, idx_z,:]
end

# ╔═╡ 1e011fdb-1d6c-471d-9d81-630be46de1c8
begin
	sens_julia = espirit(crop(kbart,(24,24,24)),(128,128,128),(6,6,6),nmaps=1,eigThresh_2=0.0);
	
	t_julia = @elapsed sens_julia = espirit(crop(kbart,(24,24,24)),(128,128,128),(6,6,6),nmaps=1,eigThresh_2=0.0);
end

# ╔═╡ cc6de5aa-d4cb-4f73-b648-1c6b8969cf98
md"# Benchmark results"

# ╔═╡ fb569c34-49d7-4954-b559-d65dc4331c82
begin
	t_bart2 = Int(round(t_bart))
	t_julia2 = Int(round(t_julia))
end;

# ╔═╡ 7466acc1-34e3-4112-b8f4-57e078a0305e
md"### BART = $t_bart2 s vs JULIA = $t_julia2 s"

# ╔═╡ 5877d93b-b2bd-4b73-b611-a5dc658aa35a
	begin
		p_vec = Any[]
		for i in 1:8
		push!(p_vec,heatmap(abs.(sens_bart[:,:,64,i]), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,showaxis = false))
		push!(p_vec,heatmap(abs.(sens_julia[:,:,64,i,1]), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,showaxis = false))
		end
		plot(p_vec...,layouts=(8,2),size=(500,1000))
			
	end

# ╔═╡ Cell order:
# ╠═e979c3be-3f37-11ed-13f0-d966a2fe1aba
# ╠═b9cf2100-3959-454a-9ba5-d49bbb0c826c
# ╠═b25c25f1-3154-40f5-a58d-f945a29ed5f8
# ╠═8087d5a8-99eb-4c11-97f2-883e41bce5b6
# ╠═550c0b48-0d94-4bf2-a15b-0cebcffe58bb
# ╠═3951341a-f2e7-465e-b783-a92d56bc05e1
# ╠═c7295363-31e2-4053-8d8c-abeed4274a89
# ╠═1e011fdb-1d6c-471d-9d81-630be46de1c8
# ╟─cc6de5aa-d4cb-4f73-b648-1c6b8969cf98
# ╠═fb569c34-49d7-4954-b559-d65dc4331c82
# ╟─7466acc1-34e3-4112-b8f4-57e078a0305e
# ╠═5877d93b-b2bd-4b73-b611-a5dc658aa35a
