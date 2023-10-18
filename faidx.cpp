#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <iomanip>
#include <cstdlib>

namespace fs = std::filesystem;

void get_file_paths(char const *fname, std::string& fasta_pathname, std::string& faidx_pathname)
{
    fs::path fasta_path = fname;

    if (!fs::exists(fasta_path))
    {
        std::cerr << "error: path " << std::quoted(fname) << " does not exist" << std::endl;
        exit(-1);
    }

    if (!fs::is_regular_file(fasta_path))
    {
        std::cerr << "error: path " << std::quoted(fname) << " does not reference a file" << std::endl;
        exit(-1);
    }

    fs::path faidx_path = std::string(fname) + ".fai";

    if (fs::exists(faidx_path))
    {
        std::cerr << "path " << std::quoted(faidx_path.string()) << " already exists, overwriting" << std::endl;
        fs::remove(faidx_path);
    }

    fasta_pathname = fasta_path.string();
    faidx_pathname = faidx_path.string();
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "usage: " << argv[0] << " <reads.fa>" << std::endl;
        return -1;
    }

    std::string fasta_pathname, faidx_pathname;
    get_file_paths(argv[1], fasta_pathname, faidx_pathname);

    std::ifstream fasta_stream(fasta_pathname);
    std::ofstream faidx_stream(faidx_pathname);

    std::string name, line;
    size_t fpos, foffset, readlen, linelen;
    bool expect_name;

    fpos = foffset = readlen = linelen = 0;
    expect_name = false;

    size_t linenum = 0;

    while (std::getline(fasta_stream, line))
    {
        line += "\n";
        linenum++;

        if (expect_name && line[0] != '>')
        {
            std::cerr << "Format error: unexpected " << std::quoted(std::string(1, line[0])) << " at line " << linenum << std::endl;
            fasta_stream.close();
            faidx_stream.close();
            fs::path p = faidx_pathname;
            fs::remove(p);
            exit(-1);
        }

        if (line[0] == '>')
        {
            if (readlen > 0)
            {
                faidx_stream << name << "\t" << readlen << "\t" << fpos << "\t" << linelen-1 << "\t" << linelen << "\n";
                readlen = 0;
            }

            name = line.substr(1, line.find_first_of(" \t\n\r\v", 1) - 1);
            expect_name = false;
        }
        else
        {
            if (readlen == 0)
            {
                fpos = foffset;
                linelen = line.size();
            }
            else if (line.size() > linelen)
            {
                std::cerr << "Format error: different line length in sequence " << std::quoted(name) << " at line " << linenum << std::endl;
                fasta_stream.close();
                faidx_stream.close();
                fs::path p = faidx_pathname;
                fs::remove(p);
                exit(-1);
            }
            else if (line.size() < linelen)
            {
                expect_name = true;
            }
            readlen += line.size()-1;
        }
        foffset += line.size();
    }

    if (readlen > 0)
    {
        faidx_stream << name << "\t" << readlen << "\t" << fpos << "\t" << linelen-1 << "\t" << linelen << "\n";
    }

    fasta_stream.close();
    faidx_stream.close();

    return 0;
}
