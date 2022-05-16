import sys

if __name__ == "__main__":
    filename = sys.argv[1]
    contents = "".join([line for line in open(filename, "r")])

    ptr = contents.find("Main:")
    time = float(contents[ptr:].split(" ", 1)[0].split(":")[1])/1000
    print("Main = {}".format(round(time,2)))

    ptr = contents.find("Dfd:PfrReadFasta()")
    time = float(contents[ptr:].split(" ", 1)[0].split(":")[2])/1000
    print("Main:PfrReadFasta() = {}".format(round(time,2)))

    ptr = contents.find("KmerOp:GenerateA:FirstPass()")
    time = float(contents[ptr:].split(" ", 1)[0].split(":")[3])/1000
    ptr = contents.find("KmerOp:GenerateA:SecondPass()")
    time += float(contents[ptr:].split(" ", 1)[0].split(":")[3])/1000
    print("KmerOp:GenerateA:FirstPass() + KmerOp:GenerateA:SecondPass() = {}".format(round(time,2)))

    ptr = contents.find("KmerOp:GenerateA:SpMatA()")
    time = float(contents[ptr:].split(" ", 1)[0].split(":")[3])/1000
    print("KmerOp:GenerateA:SpMatA() = {}".format(round(time,2)))

    ptr = contents.find("Main:AAt()")
    time = float(contents[ptr:].split(" ", 1)[0].split(":")[2])/1000
    print("Main:AAt() = {}".format(round(time,2)))

    ptr = contents.find("Main:DfdWait()")
    time = float(contents[ptr:].split(" ", 1)[0].split(":")[2])/1000
    print("Main:DdfWait() = {}".format(round(time,2)))

    ptr = contents.find("Main:DprAlign()")
    time = float(contents[ptr:].split(" ", 1)[0].split(":")[2])/1000
    print("Main:DprAlign() = {}".format(round(time,2)))

    ptr = contents.find("Main:TransitiveReduction()")
    time = float(contents[ptr:].split(" ", 1)[0].split(":")[2])/1000
    print("Main:TransitiveReduction() = {}".format(round(time,5)))

    ptr = contents.find("Main:ExtractContig()")
    time = float(contents[ptr:].split(" ", 1)[0].split(":")[2])/1000
    print("Main:ExtractContig() = {}".format(round(time,5)))

    ptr = contents.find("CreateContig:GetRead2Contigs()")
    time = float(contents[ptr:].split(" ", 1)[0].split(":")[2])/1000
    print("CreateContig:GetRead2Contigs() = {}".format(round(time,5)))

    ptr = contents.find("CreateContig:GetRead2ProcAssignments()")
    time = float(contents[ptr:].split(" ", 1)[0].split(":")[2])/1000
    print("CreateContig:GetRead2ProcAssignments()() = {}".format(round(time,5)))

    ptr = contents.find("CreateContig:InducedSubgraphs2Procs()")
    time = float(contents[ptr:].split(" ", 1)[0].split(":")[2])/1000
    print("CreateContig:InducedSubgraphs2Procs() = {}".format(round(time,5)))

    ptr = contents.find("CreateContig:ReadExchange()")
    time = float(contents[ptr:].split(" ", 1)[0].split(":")[2])/1000
    print("CreateContig:ReadExchange() = {}".format(round(time,5)))

    ptr = contents.find("CreateContig:LocalAssembly()")
    time = float(contents[ptr:].split(" ", 1)[0].split(":")[2])/1000
    print("CreateContig:LocalAssembly() = {}".format(round(time,5)))

    ptr = contents.find("Main:WriteContigs()")
    time = float(contents[ptr:].split(" ", 1)[0].split(":")[2])/1000
    print("Main:WriteContigs() = {}".format(round(time,5)))
    
