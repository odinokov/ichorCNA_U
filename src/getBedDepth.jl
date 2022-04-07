using QuickArgParse
using CSV,DataFrames,Query

function main()
	req = ["InputFile","BinSize","ChromSize","OutFile"]
	reqhelp = ["Path to input BED file","Bin size in kb",
	"Tab-delimited file of chrom sizes","Output wig file path"]
	opts = ["toolpath","chroms"]
	optsHelp=["Path to bedtools",
	"Chromosomes to analyze. Comma separated, can include ranges (ie 1:22,X). ALL for all in chrom size file"]
	optsDefault=["bedtools","ALL"]
	flags = ["p"]
	flagHelp = ["Add 'chr' prefix to all chrom names."]
	title = "BED Bin Read Counter"
	desc = "Counts total reads in fixed size bins from a bed file. 
Generates a temp bed file of the bins. Requires bedtools."

	R = process_reqs(req;reqHelp=reqhelp,title=title,desc=desc,
		optTag=opts,optHelp=optsHelp,optDefault=optsDefault,
		flags=flags,flagHelp=flagHelp)
	build_usage!(R)
	A = parse_args(R)

	inBed = A["InputFile"]
	outName = A["OutFile"]
	toolPath = A["toolpath"]
	if occursin(".wig",outName) == false
		outName = outName * ".wig"
	end

	cSize = readlines(open(A["ChromSize"],"r"))
	So = String[]           #Order of chroms in the size file
	S = Dict{String,Int}()  #Size lookup dict
	for ln in cSize
		if length(ln) == 0
			continue
		end
		data = [string(x) for x in split(ln,"\t")]
		push!(So,data[1])
		S[data[1]] = parse(Int,data[2])
	end

	#Process chrom size input file and generate bins
	binSize = parse(Int,A["BinSize"]) * 1000
	tempBed = replace(basename(outName),".wig"=>"") * ".ichorCNAtempFixedBins.bed"
	bedF = open(tempBed,"w")
	for c in So
		allS = collect(range(1,stop=S[c],step=binSize))
		allE = [x + binSize for x in allS]
		bins = [ [allS[x],allE[x]] for x in 1:length(allS) ]
		for b in bins
			write(bedF,"$c\t$(b[1])\t$(b[2])\n")
		end
	end
	close(bedF)

	#Process chromosome inputs
	if A["chroms"] == "ALL"
		allChrom = So
	else
		allChrom = String[]
		rawChrom = [string(x) for x in split(A["chroms"],",")]
		rF = r"\d+\:\d+"
		for c in rawChrom
			if isnothing(match(rF,c)) == false && length(split(c,":")) == 2
				cR = [parse(Int,x) for x in split(c,":")]
				allC = [string(x) for x in collect(range(cR[1],stop=cR[2]))]
				if A["p"] == true
					allC = ["chr" * x for x in allC]
				end
				for c in allC
					push!(allChrom,c)
				end
			else
				if A["p"] == true
					c = "chr" * c
					push!(allChrom,c)
				else
					push!(allChrom,c)
				end
			end
		end
	end

	cmd = `$toolPath coverage -counts -sorted -g $(A["ChromSize"]) -a $tempBed -b $inBed`
	colTypes = Dict(1=>String,4=>Int)
	countData = CSV.File(open(cmd,"r"),delim='\t',types=colTypes,header=false) |> DataFrame
	
	outF = open(outName,"w")
	for c in allChrom
		cData = countData |> @filter(_.Column1 == c) |> DataFrame
		write(outF,"fixedStep chrom=$c start=1 step=$binSize span=$binSize\n")
		for r in eachrow(cData)
			write(outF,"$(r[:Column4])\n")
		end
	end
	close(outF)
	run(`rm $tempBed`)

	return
end

main()
