#include "EnuBundle.h"
#include "args.hxx"
#define FILELEN 1024
// #pragma comment(linker, "/STACK:102400000,102400000")
//#define TRANSFER
//#include<vld.h>
void showUsage() {
	fprintf(stderr, "enplex [-f filename] [--small] [-k k] [-q lb] [-t maxsecond] \n");
}

int *twoPow;

int main(int argc, char** argv) {

	twoPow = new int[1 << 16];

	for (int i = 0; i < 16; ++i)
		twoPow[1 << i] = i;

	//p2p-Gnutella04,wiki-vote.txt
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\10th_dimacs\\jazz.bin";
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\snap\\soc-Slashdot0902.bin";
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\snap\\wiki-Vote.bin";
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\snap\\email-EuAll.bin";
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\snap\\amazon0505.bin";
	//char filepath[1024] = "D:\\Home\\benchmarks\\splex\\10th_dimacs\\celegans_metabolic.bin";
	char filepath[1024] = "D:\\Home\\benchmarks\\splex\\snap\\as-caida20071105.bin";
	//char filepath[1024] = "graph1.bin";
	//char filepath[1024] = "graph2.bin";
	//char filepath[FILELEN] = "D:\\Home\\benchmarks\\splex\\2nd_dimacs\\brock200_1.bin";

	ui k = 2;
	ui lb = 1;
	uli maxsec = 600;
	ui decompose = 0;
	ui isquiete = 0;

	args::ArgumentParser parser(
        "Enplex, a software for enumerating kplex\n");

    args::HelpFlag help(parser, "help", "Display this help menu",
                        {'h', "help"});
    args::Group required(parser, "", args::Group::Validators::All);

    args::ValueFlag<std::string> benchmark_file(
        parser, "benchmark", "Path to benchmark", {'f', "file"},
        "");

	args::ValueFlag<int> K(parser, "para k", "The parameter k", {'k', "k"}, 2);

    args::ValueFlag<int> TimeLim(parser, "Time limitation",
                                 "The cut down time in second", {'t', "time"},
                                 90);

	args::ValueFlag<int> LowerBound(parser, "Lower Bound", "The lower bound of the size of kplex", {'l', "lower"}, 1);

	args::ValueFlag<int> Decompose(parser, "decompose", "Decompose or not", {'d', "d"}, 0);

	args::ValueFlag<int> Quiete(parser, "quiete", "quiete or not", {'q', "q"}, 0);

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 0;
    } catch (args::ValidationError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 0;
    }

	strncpy(filepath, args::get(benchmark_file).c_str(), FILELEN);
	k = args::get(K);
	maxsec = args::get(TimeLim);
	lb = args::get(LowerBound);
	decompose = args::get(Decompose);
	isquiete = args::get(Quiete);

	if (decompose && lb < 2*k-2) {
		fprintf(stderr, "lb is at least 2k-2 in decompose mode\n");
		exit(-1);
	}
	EnuBundle enbundle;
	enbundle.readBinaryGraph(filepath);
	enbundle.enumPlex(k,lb,maxsec, decompose,isquiete);
	//_CrtDumpMemoryLeaks();
	return 0;
}

