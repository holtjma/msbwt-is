#include <string>

using namespace std;

struct Alphabet {
    unsigned long long numSym;
    unsigned long long numTerm;
};

void preprocessStrings(string bwtDir, vector<string> seqs);
template <class T> vector<unsigned long long> scanSeq(T * valArray, unsigned long long valArrayLen, unsigned long numSym, unsigned long numTerm);
void alternateBuilder(string bwtDir);
template <class T> Alphabet recurseUp(string bwtDir, unsigned long numSym, unsigned long numTerm, unsigned long recursionLevel);
template <class T> Alphabet induce(T * converted, unsigned long long conLen, vector<unsigned long long> &sStarArray, string bwtDir, unsigned long recursionLevel, unsigned long long numSym, unsigned long long numTerm);

template <class T, class U> void recurseDown(string bwtDir, Alphabet currAlphabet, Alphabet higherAlphabet, unsigned long recursionLevel);

template <class T, class U> void recurseDown2(string bwtDir, Alphabet currAlphabet, Alphabet higherAlphabet, unsigned long recursionLevel);

template <class T, class U> void recurseDown3(string bwtDir, Alphabet currAlphabet, Alphabet higherAlphabet, unsigned long recursionLevel);