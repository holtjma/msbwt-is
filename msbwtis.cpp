#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "msbwtis.h"

using namespace std;

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

inline void setBit64(vector<unsigned long long> &bitArray, unsigned long long index) {
    bitArray[index >> 6] |= ((unsigned long long)0x1) << (index & 0x3F);
}

inline bool getBit64(vector<unsigned long long> &bitArray, unsigned long long index) {
    return (bitArray[index >> 6] >> (index & 0x3F)) & 0x1;
}

template <class T> class SStarSubstring {
public:
    T * seq;
    T finalSym;
    unsigned long long l;
    
    bool operator< (const SStarSubstring &rhs) const {
        if(l < rhs.l) {
            for(unsigned long long i = 0; i < l; i++) {
                if(seq[i] != rhs.seq[i]) {
                    return (seq[i] < rhs.seq[i]);
                }
            }
            return (finalSym < rhs.seq[l]);
        }
        else if(l > rhs.l) {
            for(unsigned long long i = 0; i < rhs.l; i++) {
                if(seq[i] != rhs.seq[i]) {
                    return (seq[i] < rhs.seq[i]);
                }
            }
            return (seq[rhs.l] <= rhs.finalSym);
        }
        else {
            for(unsigned long long i = 0; i < l; i++) {
                if(seq[i] != rhs.seq[i]) {
                    return (seq[i] < rhs.seq[i]);
                }
            }
            return (finalSym < rhs.finalSym);
        }
    }
    
    bool operator== (const SStarSubstring &rhs) const {
        if(l == rhs.l) {
            for(unsigned long long i = 0; i < l; i++) {
                if(seq[i] != rhs.seq[i]) {
                    return false;
                }
            }
            return (finalSym == rhs.finalSym);
        }
        else {
            return false;
        }
    }
};

namespace std {
    template <class T> struct hash<SStarSubstring<T> >
    {
        size_t operator() (SStarSubstring<T> const& tup) const {
            /*
             newer one here
             http://www.cse.yorku.ca/~oz/hash.html - sdbm
             */
            ///*
            //time wise, there seems to be little difference between 1.1 and 1.2; however, I believe 1.1 is the most "correct" way to do this
            //1.1
            unsigned char * p = (unsigned char*)tup.seq;
            size_t ret = 0;
            for(unsigned long long i = 0; i < sizeof(T)*tup.l; i++) {
                ret = p[i] + (ret << 6) + (ret << 16) - ret;
            }
            p = (unsigned char*)&tup.finalSym;
            for(unsigned long long i = 0; i < sizeof(T); i++) {
                ret = p[i] + (ret << 6) + (ret << 16) - ret;
            }
            return ret;
            //*/
            /*
             //1.2
             size_t ret = 0;
             for(unsigned long long i = 0; i < tup.l; i++) {
             ret = tup.seq[i] + (ret << 6) + (ret << 16) - ret;
             }
             ret = tup.finalSym + (ret << 6) + (ret << 16) - ret;
             return ret;
             //*/
        }
    };
}

void preprocessStrings(string bwtDir, vector<string> seqs) {
    //by default, all characters are 'N', aka 4
    vector<char> charMap = vector<char>(256, 4);
    charMap['$'] = 0;
    charMap['A'] = 1;
    charMap['C'] = 2;
    charMap['G'] = 3;
    charMap['N'] = 4;
    charMap['T'] = 5;
    
    string seq;
    char val;
    
    string filename = bwtDir+"/seqs_test.0.dat";
    FILE * fp = fopen(filename.c_str(), "w");
    for(int i = 0; i < seqs.size(); i++) {
        seq = seqs[i];
        val = 0;
        fwrite(&val, 1, 1, fp);
        for(int j = 0; j < seq.length(); j++) {
            val = charMap[seq[j]];
            fwrite(&val, 1, 1, fp);
        }
    }
    
    fclose(fp);
}

void preprocessFasta(string bwtDir, vector<string> fastaFNs) {
    //by default, all characters are 'N', aka 4
    vector<char> charMap = vector<char>(256, 4);
    charMap['$'] = 0;
    charMap['A'] = 1;
    charMap['C'] = 2;
    charMap['G'] = 3;
    charMap['N'] = 4;
    charMap['T'] = 5;
    
    string line;
    char val = 0;
    bool prevLabel = false;
    
    string filename = bwtDir+"/seqs_test.0.dat";
    FILE * fp = fopen(filename.c_str(), "w");
    for(int i = 0; i < fastaFNs.size(); i++) {
        string fn = fastaFNs[i];
        ifstream ifp(fn);
        
        while(getline(ifp, line)) {
            if(line[0] == '>') {
                prevLabel = true;
            }
            else {
                if(prevLabel) {
                    fwrite(&val, 1, 1, fp);
                    prevLabel = false;
                }
                for(int x = 0; x < line.length(); x++) {
                    line[x] = charMap[line[x]];
                }
                fwrite(line.c_str(), 1, line.length(),fp);
            }
        }
        ifp.close();
    }
    fclose(fp);
}

template <class T> vector<unsigned long long> scanSeq(string bcFN, T * valArray, unsigned long long valArrayLen, unsigned long numSym, unsigned long numTerm) {
    vector<unsigned long long> ret = vector<unsigned long long>(ceil(valArrayLen/64.0), 0);
    bool cmp1 = false, cmp2 = false;
    
    vector<uint64_t> bc = vector<uint64_t>(numSym, 0);
    bc[valArray[valArrayLen-1]] += 1;
    
    for(int64_t x=valArrayLen-2; x >=0; x--) {
        cmp2 = cmp1;
        cmp1 = (valArray[x] < valArray[x+1]) || (valArray[x] == valArray[x+1] && cmp2);
        if(valArray[x+1] < numTerm || (!cmp1 && cmp2)) {
            setBit64(ret, x+1);
        }
        bc[valArray[x]] += 1;
    }
    setBit64(ret, 0);
    
    //save the bincount before returning
    ofstream bcstream(bcFN, ios::out | ofstream::binary);
    bcstream.write((char*)&bc[0], 8*bc.size());
    bcstream.close();
    
    return ret;
}

void alternateBuilder(string bwtDir) {
    vector<Alphabet> myAlphabets = vector<Alphabet>();
    Alphabet currAlphabet;
    currAlphabet.numSym = 6;
    currAlphabet.numTerm = 1;
    myAlphabets.push_back(currAlphabet);
    
    //here is the upwards recursion
    unsigned long recursionLevel = 0;
    while(currAlphabet.numTerm < currAlphabet.numSym) {
        if(currAlphabet.numSym <= 256) {
            currAlphabet = recurseUp<uint8_t>(bwtDir, currAlphabet.numSym, currAlphabet.numTerm, recursionLevel);
        }
        else if(currAlphabet.numSym <= pow(256, 2)) {
            currAlphabet = recurseUp<uint16_t>(bwtDir, currAlphabet.numSym, currAlphabet.numTerm, recursionLevel);
        }
        else if(currAlphabet.numSym <= pow(256, 4)) {
            currAlphabet = recurseUp<uint32_t>(bwtDir, currAlphabet.numSym, currAlphabet.numTerm, recursionLevel);
        }
        else {
            //???
        }
        myAlphabets.push_back(currAlphabet);
        recursionLevel++;
    }
    
    //whatever level we are at can be trivially solved with a sort
    string seqFN = bwtDir+"/seqs_test."+to_string(recursionLevel)+".dat";
    string bwtFN = bwtDir+"/bwt_test."+to_string(recursionLevel)+".dat";
    rename(seqFN.c_str(), bwtFN.c_str());
    
    struct stat st;
    stat(bwtFN.c_str(), &st);
    int fd = open(bwtFN.c_str(), O_RDWR, 0);
    
    //load the info from disk
    if(currAlphabet.numSym <= 256) {
        uint8_t * converted = (uint8_t *)mmap(NULL, st.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        sort(converted, converted+st.st_size);
        munmap(converted, st.st_size);
    }
    else if(currAlphabet.numSym <= pow(256, 2)) {
        uint16_t * converted = (uint16_t *)mmap(NULL, st.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        sort(converted, converted+st.st_size/2);
        munmap(converted, st.st_size);
    }
    else if(currAlphabet.numSym <= pow(256, 4)) {
        uint32_t * converted = (uint32_t *)mmap(NULL, st.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        sort(converted, converted+st.st_size/4);
        munmap(converted, st.st_size);
    }
    else {
        //TODO: currently, we don't handle alphabet this big
    }
    close(fd);
    
    while(recursionLevel > 0) {
        recursionLevel--;
        
        if(myAlphabets[recursionLevel].numSym <= 256) {
            if(myAlphabets[recursionLevel+1].numSym <= 256) {
                recurseDown<uint8_t, uint8_t>(bwtDir, myAlphabets[recursionLevel], myAlphabets[recursionLevel+1], recursionLevel);
            }
            else if(myAlphabets[recursionLevel+1].numSym <= pow(256, 2)) {
                recurseDown<uint8_t, uint16_t>(bwtDir, myAlphabets[recursionLevel], myAlphabets[recursionLevel+1], recursionLevel);
            }
            else if(myAlphabets[recursionLevel+1].numSym <= pow(256, 4)) {
                recurseDown<uint8_t, uint32_t>(bwtDir, myAlphabets[recursionLevel], myAlphabets[recursionLevel+1], recursionLevel);
            }
            else {
                printf("UNHANDLED\n");
                break;
            }
        }
        else if(myAlphabets[recursionLevel].numSym <= pow(256, 2)) {
            if(myAlphabets[recursionLevel+1].numSym <= 256) {
                recurseDown<uint16_t, uint8_t>(bwtDir, myAlphabets[recursionLevel], myAlphabets[recursionLevel+1], recursionLevel);
            }
            else if(myAlphabets[recursionLevel+1].numSym <= pow(256, 2)) {
                recurseDown<uint16_t, uint16_t>(bwtDir, myAlphabets[recursionLevel], myAlphabets[recursionLevel+1], recursionLevel);
            }
            else if(myAlphabets[recursionLevel+1].numSym <= pow(256, 4)) {
                recurseDown<uint16_t, uint32_t>(bwtDir, myAlphabets[recursionLevel], myAlphabets[recursionLevel+1], recursionLevel);
            }
            else {
                printf("UNHANDLED\n");
                break;
            }
        }
        else if(myAlphabets[recursionLevel].numSym <= pow(256, 4)) {
            if(myAlphabets[recursionLevel+1].numSym <= 256) {
                recurseDown<uint32_t, uint8_t>(bwtDir, myAlphabets[recursionLevel], myAlphabets[recursionLevel+1], recursionLevel);
            }
            else if(myAlphabets[recursionLevel+1].numSym <= pow(256, 2)) {
                recurseDown<uint32_t, uint16_t>(bwtDir, myAlphabets[recursionLevel], myAlphabets[recursionLevel+1], recursionLevel);
            }
            else if(myAlphabets[recursionLevel+1].numSym <= pow(256, 4)) {
                recurseDown<uint32_t, uint32_t>(bwtDir, myAlphabets[recursionLevel], myAlphabets[recursionLevel+1], recursionLevel);
            }
            else {
                printf("UNHANDLED\n");
                break;
            }
        }
        else {
            printf("UNHANDLED\n");
            break;
        }
        
        printf("Finished recursion at level %ld\n\n", recursionLevel);
    }
}

template <class T> Alphabet recurseUp(string bwtDir, unsigned long numSym, unsigned long numTerm, unsigned long recursionLevel) {
    printf("Performing recursive step at level %lu\n", recursionLevel);
    printf("Input alphabet size: %lu (%lu)\n", numSym, numTerm);
    
    //get the file size
    struct stat st;
    string filename = bwtDir+"/seqs_test."+to_string(recursionLevel)+".dat";
    //string filename = bwtDir+"/seqs.1.dat";
    stat(filename.c_str(), &st);
    int fd = open(filename.c_str(), O_RDONLY, 0);
    
    //load the info from disk
    T * converted = (T *)mmap(NULL, st.st_size, PROT_READ, MAP_SHARED, fd, 0);
    
    int tsize = sizeof(T);
    
    //calculate the size in terms of number of values
    //printf("%lld, %d\n", st.st_size, tsize);
    unsigned long conLen = st.st_size / tsize;
    printf("Input collection length: %lu\n", conLen);
    
    //get the S* array
    string bcFN = bwtDir+"/bc_test."+to_string(recursionLevel)+".dat";
    vector<unsigned long long> sStarArray = scanSeq<T>(bcFN, converted, conLen, numSym, numTerm);
    
    //induce using the S* array
    Alphabet ret = induce(converted, conLen, sStarArray, bwtDir, recursionLevel, numSym, numTerm);
    
    //release the map
    munmap(converted, st.st_size);
    close(fd);
    
    printf("\n");
    
    return ret;
}

template <class T, class U> void saveNextSeq(T * converted, unordered_map<SStarSubstring<T>, unsigned long long> &myUmap, ofstream &seqOut, vector<unsigned long long> &sStarArray, unsigned long long conLen, unsigned long long numTerm) {
    
    U writeVal;
    
    SStarSubstring<T> tup;
    unsigned long prevTerm = converted[0];
    unsigned long long prev1 = 0;
    
    unsigned long sU = sizeof(U);
    
    //go through the whole collection pulling out all S* substrings and adding them to the hash
    for(unsigned long long x = 1; x < conLen; x++) {
        if(getBit64(sStarArray, x)) {
            if(converted[prev1] < numTerm) {
                prevTerm = converted[prev1];
            }
            
            if(x - prev1 > 1) {
                tup.seq = &(converted[prev1]);
                tup.l = x - prev1;
                if(converted[x] < numTerm) {
                    tup.finalSym = prevTerm;
                }
                else {
                    tup.finalSym = converted[x];
                }
                
                writeVal = myUmap[tup];
                seqOut.write((char*)&writeVal, sU);
            }
            prev1 = x;
        }
    }
    
    //add the final S* substring
    if(converted[prev1] < numTerm) {
        prevTerm = converted[prev1];
    }
    if(conLen - prev1 > 1) {
        tup.seq = &(converted[prev1]);
        tup.l = conLen - prev1;
        tup.finalSym = prevTerm;
        writeVal = myUmap[tup];
        seqOut.write((char*)&writeVal, sU);
    }
    
}

template <class T> Alphabet induce(T * converted, unsigned long long conLen, vector<unsigned long long> &sStarArray, string bwtDir, unsigned long recursionLevel, unsigned long long numSym, unsigned long long numTerm) {
    SStarSubstring<T> tup;
    
    unordered_map<SStarSubstring<T>, unsigned long long> myUmap = unordered_map<SStarSubstring<T>, unsigned long long>();
    
    unsigned long prevTerm = converted[0];
    unsigned long long prev1 = 0;
    
    double st = get_wall_time();
    
    unsigned long long nextSize = 0;
    vector<unsigned long long> looperCounts = vector<unsigned long long>(numTerm);
    
    //go through the whole collection pulling out all S* substrings and adding them to the hash
    for(unsigned long long x = 1; x < conLen; x++) {
        if(getBit64(sStarArray, x)) {
            if(converted[prev1] < numTerm) {
                prevTerm = converted[prev1];
            }
            
            if(x - prev1 == 1) {
                looperCounts[converted[prev1]] += 1;
            }
            else {
                tup.seq = &(converted[prev1]);
                tup.l = x - prev1;
                if(converted[x] < numTerm) {
                    tup.finalSym = prevTerm;
                }
                else {
                    tup.finalSym = converted[x];
                }
                
                myUmap[tup] = 0;
                nextSize++;
            }
            prev1 = x;
        }
    }
    
    //add the final S* substring
    if(converted[prev1] < numTerm) {
        prevTerm = converted[prev1];
    }
    if(conLen - prev1 == 1) {
        looperCounts[converted[prev1]] += 1;
    }
    else {
        tup.seq = &(converted[prev1]);
        tup.l = conLen - prev1;
        tup.finalSym = prevTerm;
        myUmap[tup] = 0;
        nextSize++;
    }
    
    double et = get_wall_time();
    printf("library solved in %lf\n", et-st);
    
    //pull out all the S* substrings and sort them
    vector<SStarSubstring<T> > v = vector<SStarSubstring<T> >(myUmap.size());
    unsigned long long alphabetLen = 0;
    unsigned long long i = 0;
    for(typename unordered_map<SStarSubstring<T>, unsigned long long>::iterator it = myUmap.begin(); it != myUmap.end(); it++) {
        v[i] = it->first;
        i++;
        alphabetLen += it->first.l+1;
    }
    sort(v.begin(), v.end());
    printf("sort finished in %lf\n", get_wall_time()-et);
    
    //assign their values in the hash based on their lexicographic order
    ofstream alphabetOut(bwtDir+"/alphabet_test."+to_string(recursionLevel)+".dat", ios::out | ios::binary);
    alphabetOut.seekp(sizeof(T)*alphabetLen-1, ios::beg);
    alphabetOut.write("\0", 1);
    alphabetOut.flush();
    alphabetOut.seekp(0, ios::beg);
    
    ofstream offsetsOut(bwtDir+"/offsets_test."+to_string(recursionLevel)+".dat", ios::out | ios::binary);
    offsetsOut.seekp(8*myUmap.size()-1, ios::beg);
    offsetsOut.write("\0", 1);
    offsetsOut.flush();
    offsetsOut.seekp(0, ios::beg);
    
    uint64_t offset = 0;
    offsetsOut.write((char*)&offset, sizeof(uint64_t));
    unsigned long long symCount = 0;
    unsigned long long symTerm = 0;
    for(typename vector<SStarSubstring<T> >::iterator it = v.begin(); it != v.end(); it++) {
        //set the map and count this as a new symbol
        myUmap[*it] = symCount;
        symCount++;
        if(it->seq[0] < numTerm) {
            symTerm++;
        }
        
        //save the alphabet information
        alphabetOut.write((char*)it->seq, sizeof(T)*it->l);
        alphabetOut.write((char*)&(it->finalSym), sizeof(T));
        offset += it->l+1;
        offsetsOut.write((char*)&offset, sizeof(uint64_t));
    }
    
    alphabetOut.close();
    offsetsOut.close();
    
    double et2 = get_wall_time();
    printf("alpha-proc solved in %lf\n", et2-et);
    
    ofstream looperOut(bwtDir+"/lc_test."+to_string(recursionLevel)+".dat", ios::out | ios::binary);
    for(unsigned long long x = 0; x < numTerm; x++) {
        looperOut.write((char*)&looperCounts[x], 8);
    }
    looperOut.close();
    
    ofstream seqOut(bwtDir+"/seqs_test."+to_string(recursionLevel+1)+".dat", ios::out | ios::binary);
    //out.write((char*)converted, sizeof(T));
    //printf("%d\n", converted[0]);
    
    //save next super string
    if(symCount <= 256) {
        seqOut.seekp(nextSize-1, ios::beg);
        seqOut.write("\0", 1);
        seqOut.flush();
        seqOut.seekp(0, ios::beg);
        saveNextSeq<T, uint8_t>(converted, myUmap, seqOut, sStarArray, conLen, numTerm);
    }
    else if(symCount <= pow(256, 2)) {
        seqOut.seekp(2*nextSize-1, ios::beg);
        seqOut.write("\0", 1);
        seqOut.flush();
        seqOut.seekp(0, ios::beg);
        saveNextSeq<T, uint16_t>(converted, myUmap, seqOut, sStarArray, conLen, numTerm);
    }
    else if(symCount <= pow(256, 4)) {
        seqOut.seekp(4*nextSize-1, ios::beg);
        seqOut.write("\0", 1);
        seqOut.flush();
        seqOut.seekp(0, ios::beg);
        saveNextSeq<T, uint32_t>(converted, myUmap, seqOut, sStarArray, conLen, numTerm);
    }
    else {
        //TODO: raise an exception or something? we are not currently allowing alphabets that large
    }
    
    seqOut.close();
    
    double et3 = get_wall_time();
    
    //printf("next size: %llu\n", nextSize);
    //printf("umap size: %lu\n", myUmap.size());
    printf("SS solved in %lf\n", et3-et2);
    
    Alphabet ret;
    ret.numSym = symCount;
    ret.numTerm = symTerm;
    return ret;
}

#define BREAK_POINT_MARKER 0xffffffffffffffff

template <class T, class U> void recurseDown(string bwtDir, Alphabet currAlphabet, Alphabet higherAlphabet, unsigned long recursionLevel) {
    
    printf("CALLED recurseDown\n");
    
    //load the alphabet, class T
    struct stat al_st;
    string alFN = bwtDir+"/alphabet_test."+to_string(recursionLevel)+".dat";
    stat(alFN.c_str(), &al_st);
    uint64_t alphabetLen = al_st.st_size / sizeof(T);
    ifstream alIn(alFN, ios::in | ios::binary);
    vector<T> alphabet = vector<T>(alphabetLen);
    alIn.read((char*)&alphabet[0], al_st.st_size);
    alIn.close();
    
    //load the offsets, uint64_t
    struct stat of_st;
    string ofFN = bwtDir+"/offsets_test."+to_string(recursionLevel)+".dat";
    stat(ofFN.c_str(), &of_st);
    ifstream ofIn(ofFN, ios::in | ios::binary);
    vector<uint64_t> ssOffsets = vector<uint64_t>(of_st.st_size/8);
    ofIn.read((char*)&ssOffsets[0], of_st.st_size);
    ofIn.close();
    
    //load the recursed BWT result, class U
    struct stat rb_st;
    string rbFN = bwtDir+"/bwt_test."+to_string(recursionLevel+1)+".dat";
    stat(rbFN.c_str(), &rb_st);
    uint64_t rbLen = rb_st.st_size / sizeof(U);
    //int rb_fd = open(rbFN.c_str(), O_RDONLY, 0);
    //U * recurseBWT = (U *)mmap(NULL, rb_st.st_size, PROT_READ, MAP_PRIVATE, rb_fd, 0);
    ifstream rbIn(rbFN, ios::in | ios::binary);
    vector<U> recurseBWT = vector<U>(rbLen);
    rbIn.read((char*)&recurseBWT[0], rb_st.st_size);
    rbIn.close();
    
    //memmap the looper counts, uint64_t
    struct stat lc_st;
    string lcFN = bwtDir+"/lc_test."+to_string(recursionLevel)+".dat";
    stat(lcFN.c_str(), &lc_st);
    int lc_fd = open(lcFN.c_str(), O_RDONLY, 0);
    uint64_t * lc = (uint64_t *)mmap(NULL, lc_st.st_size, PROT_READ, MAP_PRIVATE, lc_fd, 0);
    
    //memmap the bc stats, uint64_t
    struct stat bc_st;
    string bcFN = bwtDir+"/bc_test."+to_string(recursionLevel)+".dat";
    stat(bcFN.c_str(), &bc_st);
    int bc_fd = open(bcFN.c_str(), O_RDONLY, 0);
    uint64_t * bc = (uint64_t *)mmap(NULL, bc_st.st_size, PROT_READ, MAP_PRIVATE, bc_fd, 0);
    
    //figure out the input size
    struct stat seq_st;
    string seqFN = bwtDir+"/seqs_test."+to_string(recursionLevel)+".dat";
    stat(seqFN.c_str(), &seq_st);
    uint64_t totLen = seq_st.st_size / sizeof(T);
    
    vector<T> bwt = vector<T>(totLen, 0);
    
    //calculate recurseBC
    vector<uint64_t> recurseBC = vector<uint64_t>(higherAlphabet.numSym, 0);
    for(uint64_t x = 0; x < rbLen; x++) {
        recurseBC[recurseBWT[x]]++;
    }
    
    double st = get_wall_time();
    
    uint64_t currOffsetInd = 0;
    vector<uint64_t> offsetsValues = vector<uint64_t>(1, 0);
    vector<uint64_t> offsetPointers = vector<uint64_t>(alphabetLen, 0);
    
    int64_t i;
    
    //double et1 = get_wall_time();
    
    //prepare to induce sort
    vector<vector<pair<uint64_t, uint64_t> > *> ss = vector<vector<pair<uint64_t, uint64_t> > *>(currAlphabet.numSym);
    vector<vector<pair<uint64_t, uint64_t> > *> l = vector<vector<pair<uint64_t, uint64_t> > *>(currAlphabet.numSym);
    vector<vector<pair<uint64_t, uint64_t> > *> ls = vector<vector<pair<uint64_t, uint64_t> > *>(currAlphabet.numSym);
    vector<vector<pair<uint64_t, uint64_t> > *> s = vector<vector<pair<uint64_t, uint64_t> > *>(currAlphabet.numSym);
    
    for(uint64_t x = 0; x < currAlphabet.numSym; x++) {
        ss[x] = NULL;
        l[x] = NULL;
        ls[x] = NULL;
        s[x] = NULL;
    }
    
    pair<uint64_t, uint64_t> myPair;
    uint64_t ind;
    uint64_t sym;
    for(uint64_t x = 0; x < higherAlphabet.numSym; x++) {
        ind = ssOffsets[x+1]-1;
        sym = alphabet[ind];
        if(ss[sym] == NULL) {
            ss[sym] = new vector<pair<uint64_t, uint64_t> >();
        }
        myPair.first = x;
        myPair.second = ind;
        ss[sym]->push_back(myPair);
        //printf("ss[%lld] push: %lld %lld\n", sym, myPair.first, myPair.second);
    }
    
    //forward traversal of trie
    uint64_t currStart = 0;
    uint64_t currOffset;
    bool lsTrigger = false;
    
    uint64_t tempSym;
    
    uint64_t y;
    
    set<uint64_t> changedList = set<uint64_t>();
    
    for(uint64_t x = 0; x < currAlphabet.numSym; x++) {
        currOffset = currStart;
        currStart += bc[x];
        
        //assert(lsTrigger == false);
        
        //empty the L queues
        if(l[x] != NULL) {
            currOffsetInd++;
            offsetsValues.push_back(currOffset);
            
            changedList.clear();
            y = 0;
            while(y < l[x]->size()) {
                myPair = l[x][0][y];
                
                if(myPair.first != BREAK_POINT_MARKER) {
                    //this is not a break point, so it belongs to the same group
                    currOffset += recurseBC[myPair.first];
                    offsetPointers[myPair.second] = currOffsetInd;
                    
                    //compare the predecessor to pick the right queue
                    tempSym = alphabet[myPair.second-1];
                    if(tempSym < x) {
                        //push the one we retrieve to the LS queue
                        if(ls[x] == NULL) {
                            ls[x] = new vector<pair<uint64_t, uint64_t> >();
                        }
                        ls[x]->push_back(myPair);
                        lsTrigger = true;
                    }
                    else {
                        //push the one preceding to the LS queue
                        myPair.second--;
                        if(l[tempSym] == NULL) {
                            l[tempSym] = new vector<pair<uint64_t, uint64_t> >();
                        }
                        changedList.insert(tempSym);
                        l[tempSym]->push_back(myPair);
                    }
                }
                else {
                    //TODO: this might be making more than necessary for the very last value set from each queue
                    //push a new offset
                    currOffsetInd++;
                    offsetsValues.push_back(currOffset);
                    
                    //write our new break point markers
                    myPair.first = BREAK_POINT_MARKER;
                    for(set<uint64_t>::iterator setit = changedList.begin(); setit != changedList.end(); setit++) {
                        l[*setit]->push_back(myPair);
                    }
                    changedList.clear();
                    if(lsTrigger) {
                        ls[x]->push_back(myPair);
                        lsTrigger = false;
                    }
                }
                
                y++;
            }
            
            delete l[x];
            l[x] = NULL;
        }
        
        //empty the S* queues
        if(ss[x] != NULL) {
            changedList.clear();
            y = 0;
            while(y < ss[x]->size()) {
                myPair = ss[x][0][y];
                myPair.second--;
                tempSym = alphabet[myPair.second];
                
                if(l[tempSym] == NULL) {
                    l[tempSym] = new vector<pair<uint64_t, uint64_t> >();
                }
                changedList.insert(tempSym);
                l[tempSym]->push_back(myPair);
                y++;
            }
            
            myPair.first = BREAK_POINT_MARKER;
            for(set<uint64_t>::iterator setit = changedList.begin(); setit != changedList.end(); setit++) {
                l[*setit]->push_back(myPair);
            }
            delete ss[x];
            ss[x] = NULL;
        }
    }
    
    double et2 = get_wall_time();
    printf("Finished forward traversal in %lf seconds\n", (et2-st));
    
    //add the end of the above loop, currOffsetInd is the index of the last value in the vector, so we should increment it here
    currOffsetInd++;
    
    //reverse induce traversal
    for(i = currAlphabet.numSym-1; i >= 0; i--) {
        currOffset = currStart;
        currStart -= bc[i];
        
        //empty the S* queue
        if(s[i] != NULL) {
            changedList.clear();
            y = 0;
            while(y < s[i]->size()) {
                myPair = s[i][0][y];
                
                if(myPair.first != BREAK_POINT_MARKER) {
                    currOffset -= recurseBC[myPair.first];
                    offsetPointers[myPair.second] = currOffsetInd;
                    
                    if(ssOffsets[myPair.first] != myPair.second) {
                        myPair.second--;
                        tempSym = alphabet[myPair.second];
                        if(s[tempSym] == NULL) {
                            s[tempSym] = new vector<pair<uint64_t, uint64_t> >();
                        }
                        changedList.insert(tempSym);
                        s[tempSym]->push_back(myPair);
                    }
                    else {
                        //do nothing
                    }
                }
                else {
                    offsetsValues.push_back(currOffset);
                    currOffsetInd++;
                    
                    //push a break point now
                    myPair.first = BREAK_POINT_MARKER;
                    for(set<uint64_t>::iterator setit = changedList.begin(); setit != changedList.end(); setit++) {
                        s[*setit]->push_back(myPair);
                    }
                    changedList.clear();
                }
                
                y++;
            }
            
            delete s[i];
            s[i] = NULL;
        }
        
        offsetsValues.push_back(currOffset);
        currOffsetInd++;
        
        //empty the LS queue
        if(ls[i] != NULL) {
            changedList.clear();
            y = ls[i]->size();
            while(y > 0) {
                y--;
                myPair = ls[i][0][y];
                
                if(myPair.first != BREAK_POINT_MARKER) {
                    myPair.second--;
                    tempSym = alphabet[myPair.second];
                    if(s[tempSym] == NULL) {
                        s[tempSym] = new vector<pair<uint64_t, uint64_t> >();
                    }
                    changedList.insert(tempSym);
                    s[tempSym]->push_back(myPair);
                }
                else {
                    //write our new break point markers
                    myPair.first = BREAK_POINT_MARKER;
                    for(set<uint64_t>::iterator setit = changedList.begin(); setit != changedList.end(); setit++) {
                        s[*setit]->push_back(myPair);
                    }
                    changedList.clear();
                }
            }
            
            //do this one last time
            myPair.first = BREAK_POINT_MARKER;
            for(set<uint64_t>::iterator setit = changedList.begin(); setit != changedList.end(); setit++) {
                s[*setit]->push_back(myPair);
            }
            changedList.clear();
            
            delete ls[i];
            ls[i] = NULL;
        }
        
        //write the looper counts now
        if(i < currAlphabet.numTerm) {
            for(uint64_t x = 0; x < lc[i]; x++) {
                currOffset--;
                bwt[currOffset] = i;
            }
        }
    }
    
    //clear these out
    ss = vector<vector<pair<uint64_t, uint64_t> > *>();
    l = vector<vector<pair<uint64_t, uint64_t> > *>();
    ls = vector<vector<pair<uint64_t, uint64_t> > *>();
    s = vector<vector<pair<uint64_t, uint64_t> > *>();
    
    double et3 = get_wall_time();
    printf("Finished reverse traversal in %lf seconds\n", (et3-et2));
    
    //recurse offset computation
    vector<uint64_t> recurseOffsets = vector<uint64_t>(higherAlphabet.numSym+1, 0);
    currOffset = 0;
    for(uint64_t x = 0; x < higherAlphabet.numSym; x++) {
        recurseOffsets[x] = currOffset;
        currOffset += recurseBC[x];
    }
    recurseOffsets[higherAlphabet.numSym] = currOffset;
    
    double et4 = get_wall_time();
    printf("Finished recurse offset computation in %lf seconds\n", (et4-et3));
    
    /*
     //uncomment for debugging the output
     printf("recurseBWT\n");
     for(uint64_t ind = 0; ind < rbLen; ind++) {
     printf("%llu ", recurseBWT[ind]);
     }
     printf("\n");
     
     printf("recurseOffsets\n");
     for(uint64_t ind = 0; ind < higherAlphabet.numSym+1; ind++) {
     printf("%llu ", recurseOffsets[ind]);
     }
     printf("\n");
     
     printf("joinedSSList\n");
     for(uint64_t ind = 0; ind < alphabetLen; ind++) {
     printf("%llu ", alphabet[ind]);
     }
     printf("\n");
     
     printf("offsetPointers\n");
     for(uint64_t ind = 0; ind < offsetPointers.size(); ind++) {
     printf("%llu ", offsetPointers[ind]);
     }
     printf("\n");
     */
    
    //BWT assignment loop
    uint64_t currSym;
    uint64_t predSym;
    uint64_t finalPred;
    uint64_t runLen;
    uint64_t offTemp, offTemp2;
    
    y = 0;
    while(y < rbLen) {
        currSym = recurseBWT[y];
        runLen = 1;
        y++;
        while(y < rbLen && recurseBWT[y] == currSym) {
            runLen++;
            y++;
        }
        
        for(int64_t i = ssOffsets[currSym+1]-2; i > ssOffsets[currSym]; i--) {
            //set all values in bwt[offsetsValues[offsetPointers[i]]:+runlen) to the character before
            fill_n(bwt.begin()+offsetsValues[offsetPointers[i]], runLen, alphabet[i-1]);
            offsetsValues[offsetPointers[i]] += runLen;
        }
        
        //have to do these characters one at a time, because they may not all have the same predecessor
        offTemp = offsetsValues[offsetPointers[ssOffsets[currSym]]];
        offTemp2 = recurseOffsets[currSym];
        for(uint64_t z = 0; z < runLen; z++) {
            predSym = recurseBWT[offTemp2+z];
            finalPred = alphabet[ssOffsets[predSym+1]-2];
            bwt[offTemp+z] = finalPred;
        }
        recurseOffsets[currSym] += runLen;
        offsetsValues[offsetPointers[ssOffsets[currSym]]] += runLen;
        
        //TODO: is there any way to turn the above into a run based thing?
        //issue with this method is when one S* substring is the suffix of another; example: ATA and AATA
        //fill_n(bwt.begin()+offsetsValues[offsetPointers[ssOffsets[x]]], runLen, alphabet[ssOffsets[currSym+1]-2]);
        //offsetsValues[offsetPointers[ssOffsets[x]]] += runLen;
    }
    
    string bwtFN = bwtDir+"/bwt_test."+to_string(recursionLevel)+".dat";
    ofstream bwtOut(bwtFN, ios::out | ios::binary);
    bwtOut.write((char*)&bwt[0], sizeof(T)*totLen);
    bwtOut.close();
    
    double et5 = get_wall_time();
    printf("Finished bwt assignment in %lf seconds\n", (et5-et4));
    
    //release the various maps we have open
    remove(alFN.c_str());
    
    remove(ofFN.c_str());
    
    remove(rbFN.c_str());
    
    munmap(lc, lc_st.st_size);
    close(lc_fd);
    remove(lcFN.c_str());
    
    munmap(bc, bc_st.st_size);
    close(bc_fd);
    remove(bcFN.c_str());
    
    remove(seqFN.c_str());
    
    //don't delete this one yet
    //munmap(bwt, seq_st.st_size);
    //close(bwt_fd);
}

int main (int argc, char * argv[]) {
    //string bwtDir = "/Users/Matt/testBWT/testBWT3";
    //string bwtDir = "/Users/Matt/testBWT/LinearBuilderTester";
    
    if(argc < 3) {
        printf("Usage: ./msbwtis outputDirectory input1.fa ...\n");
        return 0;
    }
    
    string bwtDir = argv[1];
    vector<string> fastaFNs = vector<string>();
    
    //printf("%d\n", argc);
    for(int x = 2; x < argc; x++) {
        fastaFNs.push_back(string(argv[x]));
    }
    
    //fastaFNs.push_back("/Users/Matt/Downloads/pacbio.fasta");
    preprocessFasta(bwtDir, fastaFNs);
    
    //preprocessStrings(bwtDir, myList);
    alternateBuilder(bwtDir);
}