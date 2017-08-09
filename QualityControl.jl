module QualityControl

using DataFrames
using BioSequences
using BioSequences.FASTQ

export QC

const bc=["ACAG","AGTC","ATCA","CATG","CTAC","TCTA","TGAT","TTGG"]

"""
    QC(file)

counts the number of As and the longest streak of As in a FASTQ file.  Additionally, demultiplexes the reads per pool.

"""

function QC(file::String; barcodes::Vector{String}=bc)
    reader = open(FASTQ.Reader,file)

    # generate an empty DataFrame
    data = DataFrame(
        countsA=Int[],
        streakA=Int[],
        pool=String[])
    
    dplxr = Demultiplexer(
        # convert the barcodes to a proper DNA sequece
        BioSequence{DNAAlphabet{4}}[barcodes...], 
        n_max_errors=1,
        distance=:hamming)
    
    # generate a demultiplexer

    for record in reader
        seq = sequence(record)
        comp = composition(seq)
        poolid, errors = demultiplex(dplxr,seq[4:7])
        if poolid == 0
            poolname="Unknown"
        else
            poolname=barcodes[poolid]
        end
        push!(data[:countsA],comp[DNA_A])
        push!(data[:streakA],longeststreak(seq))
        push!(data[:pool],poolname)
    end

    close(reader)
    return data
end


"""

    longeststreak(seq,nuc=DNA_A)

Finds the longest streak of `nuc` in a given sequence `seq`.

"""
function longeststreak(seq; nuc=DNA_A)
    longest=0
    current=0
    for n in seq
        if n==nuc
            current+=1
            longest=max(longest,current)
        else
            current=0
        end
    end
    return longest
end

end