/**
>HEADER
    Copyright (c) 2013 -- 2017 Rob Patro rob.patro@cs.stonybrook.edu

    This file is part of salmon.

    salmon is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    salmon is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with salmon.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/
#include <backtrace.hpp>
#include <boost/thread/thread.hpp>

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/range/irange.hpp>

// C++ string formatting library
#include "spdlog/fmt/fmt.h"

#include "GenomicFeature.hpp"
#include "SalmonConfig.hpp"
#include "VersionChecker.hpp"
#include "SalmonIndex.hpp"

int help(const std::vector<std::string>& /*opts*/) {
  fmt::MemoryWriter helpMsg;
  helpMsg.write("salmon v{}\n\n", salmon::version);
  helpMsg.write(
      "Usage:  salmon -h|--help or \n"
      "        salmon -v|--version or \n"
      "        salmon -c|--cite or \n"
      "        salmon [--no-version-check] <COMMAND> [-h | options]\n\n");
  helpMsg.write("Commands:\n");
  helpMsg.write("     index      : create a salmon index\n");
  helpMsg.write("     quant      : quantify a sample\n");
  helpMsg.write("     alevin     : single cell analysis\n");
  helpMsg.write("     swim       : perform super-secret operation\n");
  helpMsg.write("     quantmerge : merge multiple quantifications into a single file\n");

  std::cout << helpMsg.str();
  return 0;
}

int dualModeMessage() {
  auto helpmsg = R"(
    ===============

    salmon quant has two modes --- one quantifies expression using raw reads
    and the other makes use of already-aligned reads (in BAM/SAM format).
    Which algorithm is used depends on the arguments passed to salmon quant.
    If you provide salmon with alignments '-a [ --alignments ]' then the
    alignment-based algorithm will be used, otherwise the algorithm for
    quantifying from raw reads will be used.

    to view the help for salmon's selective-alignment-based mode, use the command

    salmon quant --help-reads

    To view the help for salmon's alignment-based mode, use the command

    salmon quant --help-alignment

    )";
  std::cout << "    salmon v" << salmon::version << helpmsg << "\n";
  return 0;
}

typedef std::function<int(int, const char*[], std::unique_ptr<SalmonIndex>& index)> SubCmdType;

/**
 * Bonus!
 */
int salmonSwim(int /*argc*/, const char* /*argv*/[], std::unique_ptr<SalmonIndex>& /*index*/) {

  std::cout << R"(
    _____       __
   / ___/____ _/ /___ ___  ____  ____
   \__ \/ __ `/ / __ `__ \/ __ \/ __ \
  ___/ / /_/ / / / / / / / /_/ / / / /
 /____/\__,_/_/_/ /_/ /_/\____/_/ /_/


)";

  return 0;
}

/**
 * Citation
 */
void printCite() {

  std::cout << R"(
If you use salmon in your research, please cite the publication in any
papers, pre-prints or reports.  The proper citation information for salmon
appears below.

Reference:
==========

Rob Patro, Geet Duggal, Michael I. Love, Rafael A. Irizarry, Carl Kingsford.
Salmon provides fast and bias-aware quantification of transcript expression.
Nature Methods. 2017;14(4):417-419. doi: 10.1038/nmeth.4197

bibtex:
=======

@article{Patro2017Salmon,
  doi = {10.1038/nmeth.4197},
  url = {https://doi.org/10.1038%2Fnmeth.4197},
  year  = {2017},
  month = {mar},
  publisher = {{Springer Nature}},
  volume = {14},
  number = {4},
  pages = {417--419},
  author = {Rob Patro and Geet Duggal and Michael I Love and Rafael A Irizarry and Carl Kingsford},
  title = {Salmon provides fast and bias-aware quantification of transcript expression},
  journal = {{Nature Methods}}
}
)";
}

int salmonIndex(int argc, const char* argv[], std::unique_ptr<SalmonIndex>& index);
int salmonQuantify(int argc, const char* argv[], std::unique_ptr<SalmonIndex>& index);
int salmonAlignmentQuantify(int argc, const char* argv[], std::unique_ptr<SalmonIndex>& index);
int salmonAlignmentDualMode(int argc, const char* argv[], std::unique_ptr<SalmonIndex>& index);
// TODO : PF_INTEGRATION
int salmonBarcoding(int argc, const char* argv[], std::unique_ptr<SalmonIndex>& index);
int salmonQuantMerge(int argc, const char* argv[],
                     std::unique_ptr<SalmonIndex>& index);

std::vector<std::string> getStrandedness(
  const boost::filesystem::path& indexDirectory,
  std::vector<std::string> mates1_paths,
  std::vector<std::string> mates2_paths,
  std::vector<std::string> unmated_paths
);


bool verbose = false;

int main(int argc, const char* argv[]) {
	std::vector<std::string> paths_1;
	std::vector<std::string> paths_2;
	paths_1.push_back("/home/admin/salmon/sample_data/reads_1.fastq");
	paths_2.push_back("/home/admin/salmon/sample_data/reads_2.fastq");
	std::string index_path = "/home/admin/salmon_index";
	std::vector<std::string> unstranded_paths;
	// std::cout <<
	for(auto libtype :   getStrandedness(index_path, paths_1, paths_2, unstranded_paths))
		std::cout << "LibType: " << libtype  << "\n";
	return 0;
}

int salmonAlignmentDualMode(int argc, const char* argv[], std::unique_ptr<SalmonIndex>& index) {
  // If the command is quant; determine whether
  // we're quantifying with raw sequences or alignments
  if (argc < 2) {
    return dualModeMessage();
  }
  // detect mode-specific help request
  if (strncmp(argv[1], "--help-alignment", 16) == 0) {
    const char* helpArgv[] = {argv[0], "--help", nullptr};
    return salmonAlignmentQuantify(2, helpArgv, index);
  } else if (strncmp(argv[1], "--help-reads", 12) == 0) {
    const char* helpArgv[] = {argv[0], "--help", nullptr};
    return salmonQuantify(2, helpArgv, index);
  }

  // detect general help request
  if (strncmp(argv[1], "--help", 6) == 0 or strncmp(argv[1], "-h", 2) == 0) {
    return dualModeMessage();
  }

  // otherwise, detect and dispatch the correct mode
  bool useSalmonAlign{false};
  for (int i = 0; i < argc; ++i) {
    if (strncmp(argv[i], "-a", 2) == 0 or
        strncmp(argv[i], "-e", 2) == 0 or
        strncmp(argv[i], "--alignments", 12) == 0 or
        strncmp(argv[i], "--eqclasses", 11) == 0) {
      useSalmonAlign = true;
      break;
    }
  }
  if (useSalmonAlign) {
    return salmonAlignmentQuantify(argc, argv, index);
  } else {
    return salmonQuantify(argc, argv, index);
  }
}
