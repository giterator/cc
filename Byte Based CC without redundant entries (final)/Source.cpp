

#define _CRT_SECURE_NO_WARNINGS
#define ull unsigned long long int 

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>
#include<fstream>
#include <stdint.h>
#include <queue>
#include <map>
#include <climits> // for CHAR_BIT
#include <cstring>
#include <ctype.h>
#include <algorithm>
#include<bitset>
#include<string>
#include<complex>
#include<chrono>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <thread>

#include <cassert>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <ctime>

// Probabilities are expressed in fixed point, with kProbBits bits of
// resolution. No need to go overboard with this.

const int TABLE_SIZE = 256;

using namespace std;

/*

const int UniqueSymbols = 1 << CHAR_BIT;

//const int UniqueSymbols = 513;

typedef std::string HuffCode;
//typedef std::map<char, HuffCode> HuffCodeMap;
typedef std::map<int, HuffCode> HuffCodeMap;


class INode
{
public:
	const int f;

	virtual ~INode() {}

protected:
	INode(int f) : f(f) {}
};

class InternalNode : public INode
{
public:
	INode *const left;
	INode *const right;

	InternalNode(INode* c0, INode* c1) : INode(c0->f + c1->f), left(c0), right(c1) {}
	~InternalNode()
	{
		delete left;
		delete right;
	}
};

class LeafNode : public INode
{
public:
	//const unsigned char c;
	const int c;

	LeafNode(int f, unsigned char c) : INode(f), c(c) {}
};

struct NodeCmp
{
	bool operator()(const INode* lhs, const INode* rhs) const { return lhs->f > rhs->f; }
};

INode* BuildTree(const int(&frequencies)[UniqueSymbols])
{
	std::priority_queue<INode*, std::vector<INode*>, NodeCmp> trees;

	for (int i = 0; i < UniqueSymbols; ++i)
	{
		if (frequencies[i] != 0)
			trees.push(new LeafNode(frequencies[i], (unsigned char)i));
	}
	while (trees.size() > 1)
	{
		INode* childR = trees.top();
		trees.pop();

		INode* childL = trees.top();
		trees.pop();

		INode* parent = new InternalNode(childR, childL);
		trees.push(parent);
	}
	return trees.top();
}

void GenerateCodes(const INode* node, const HuffCode& prefix, HuffCodeMap& outCodes)
{
	if (const LeafNode* lf = dynamic_cast<const LeafNode*>(node))
	{
		outCodes[lf->c] = prefix;
	}
	else if (const InternalNode* in = dynamic_cast<const InternalNode*>(node))
	{
		HuffCode leftPrefix = prefix;
		//leftPrefix.push_back(false);
		leftPrefix += '0';
		GenerateCodes(in->left, leftPrefix, outCodes);

		HuffCode rightPrefix = prefix;
		//rightPrefix.push_back(true);
		rightPrefix += '1';
		GenerateCodes(in->right, rightPrefix, outCodes);
	}
}

string HuffmanEncode(vector<unsigned char> uoutput)
{
	int ufrequencies[UniqueSymbols] = { 0 };
	string ubinaryoutput;
	if (uoutput.empty() == false)
	{

		for (std::vector<unsigned char>::iterator it = uoutput.begin(); it != uoutput.end(); ++it)
		{
			++ufrequencies[static_cast<int>(*it)];
		}
		INode* uroot = BuildTree(ufrequencies);

		HuffCodeMap ucodes;
		GenerateCodes(uroot, HuffCode(), ucodes);
		delete uroot;

		for (std::vector<unsigned char>::iterator it = uoutput.begin(); it != uoutput.end(); ++it)
		{
			ubinaryoutput += ucodes.find(*it)->second;
		}
	}
	return ubinaryoutput;
	//cout << "uoutput size: " << ubinaryoutput.length() / 8 << "\n";
}
*/


class progress_bar
{
	static const auto overhead = sizeof " [100%]";

	std::ostream& os;
	const std::size_t bar_width;
	std::string message;
	const std::string full_bar;

public:
	progress_bar(std::ostream& os, std::size_t line_width,
		std::string message_, const char symbol = '.')
		: os{ os },
		bar_width{ line_width - overhead },
		message{ std::move(message_) },
		full_bar{ std::string(bar_width, symbol) + std::string(bar_width, ' ') }
	{
		if (message.size() + 1 >= bar_width || message.find('\n') != message.npos) {
			os << message << '\n';
			message.clear();
		}
		else {
			message += ' ';
		}
		write(0.0);
	}

	// not copyable
	progress_bar(const progress_bar&) = delete;
	progress_bar& operator=(const progress_bar&) = delete;

	~progress_bar()
	{
		write(1.0);
		os << '\n';
	}

	void write(double fraction);
};

void progress_bar::write(double fraction)
{
	// clamp fraction to valid range [0,1]
	if (fraction < 0)
		fraction = 0;
	else if (fraction > 1)
		fraction = 1;

	auto width = bar_width - message.size();
	auto offset = bar_width - static_cast<unsigned>(width * fraction);

	os << '\r' << message;
	os.write(full_bar.data() + offset, width);
	os << " [" << std::setw(3) << static_cast<int>(100 * fraction) << "%] " << std::flush;
}

template<std::size_t N>
void breverse(std::bitset<N>& b) {
	for (std::size_t i = 0; i < N / 2; ++i) {
		bool t = b[i];
		b[i] = b[N - i - 1];
		b[N - i - 1] = t;
	}
}

string char_to_str(unsigned char ch)
{
	bitset<8> b(ch);
	//breverse(&temp);

	/*for (std::size_t i = 0; i < 8 / 2; ++i) {
		bool t = b[i];
		b[i] = b[8 - i - 1];
		b[8 - i - 1] = t;
	}*/
	return b.to_string();
}

//concentric circles
// common letter node in hash table declaration
class HashNode
{
public:
	unsigned char c; //key is same as char
	int r1, r2, p1, p2;
	HashNode* next;
	HashNode* cnext;
	HashNode* cprev;
	HashNode(unsigned char c, int r1, int p1, int p2)
	{
		this->c = c;
		this->r1 = r1;
		this->r2 = 2; //was 0
		this->p1 = p1;
		this->p2 = p2;
		this->cnext = NULL;
		this->cprev = NULL;
		this->next = NULL;
	}
};


/*
* HashMap Class Declaration
*/
class HashMap
{
private:
	HashNode** htable;
	HashNode* head, * tail;
public:
	HashMap()
	{
		head = NULL;
		tail = NULL;
		htable = new HashNode * [TABLE_SIZE];
		for (int i = 0; i < TABLE_SIZE; i++)
			htable[i] = NULL;
	}
	~HashMap()
	{
		for (int i = 0; i < TABLE_SIZE; ++i)
		{
			HashNode* entry = htable[i];
			while (entry != NULL)
			{
				HashNode* prev = entry;
				entry = entry->next;
				delete prev;
			}
		}
		delete[] htable;
	}
	/*
	* Hash Function
	*/
	int HashFunc(unsigned char c)// digits arent chagning to their ascii decimal version so are hardcoded into conditionals
	{

		int key = static_cast<int> (c);
		return key & 255;
	}

	/*
	* Insert Element at a key at the start of the list
	*/
	void Insert(int r1, unsigned char c, int p1, int p2)
	{
		//cout << "insert : " << c << "\n";
		unsigned int hash_val = HashFunc(c);
		HashNode* prev = nullptr;
		HashNode* p = nullptr;
		HashNode* entry = nullptr;
		entry = htable[hash_val];
		prev = new HashNode(c, r1, p1, p2);
		if (head == NULL)
		{
			head = prev;
			tail = prev;
			head->cnext = NULL;
			head->cprev = NULL;

		}
		else
		{
			p = tail;
			tail = prev;
			p->cnext = tail;
			tail->cprev = p;
			tail->cnext = NULL;


		}
		htable[hash_val] = prev;
		prev->next = entry;

		//return prev;


	}

	/*void pop(HashNode *node)
	{
		HashNode *cprev = node->cprev;
		HashNode *next = node->next;
		HashNode *cnext = node->cnext;

		unsigned int hash_val = HashFunc(node->c);
		HashNode *prev = nullptr;
		//HashNode *p = nullptr;
		prev = htable[hash_val];

		if (cprev != nullptr)//node!=tail
		{
			cprev->cnext = cnext;
		}

		if (cnext != nullptr)//node!=head
		{
			cnext->cprev = cprev;
		}
		if (node == head)
		{
			head = cnext;
		}

		if (node == prev)
		{
			htable[hash_val] = next;
			//pophead(hash_val);
			//prev = prev->next;
			//htable[hash_val] = prev;
		}
		else
		{
			//prev->next = p;
			while (prev->next != node ) { prev = prev->next; }//&& prev->next != nullptr
			prev->next = next; //singly linked list
		}



		delete node;
	}*/

	void InsertBefore(int r1, unsigned char c, HashNode* index, int p1, int p2)// index wasnt here before
	{

		unsigned int hash_val = HashFunc(c);
		HashNode* prev = nullptr;
		HashNode* p = nullptr;
		HashNode* entry = nullptr;
		entry = htable[hash_val];
		prev = new HashNode(c, r1, p1, p2);

		if (index->cprev != nullptr) // prev->cprev != NULL   index!=head
		{

			p = index->cprev;


			prev->cnext = index;
			prev->cnext->cprev = prev;
			p->cnext = prev;
			prev->cprev = p;

			htable[hash_val] = prev;
			prev->next = entry;

		}
		else
		{

			p = head;
			head = prev;
			head->cnext = index;
			head->cnext->cprev = head;

			htable[hash_val] = prev;
			prev->next = entry;

		}

		//return prev;
	}

	/*void InsertAfter(int r1, unsigned char c, HashNode *index)// index wasnt here before
	{

		unsigned int hash_val = HashFunc(c);
		HashNode *prev = nullptr;
		HashNode *p = nullptr;
		HashNode *entry = nullptr;
		entry = htable[hash_val];
		prev = new HashNode(c, r1);

		p = index->cnext;


		prev->cnext = p;
		prev->cnext->cprev = prev;
		index->cnext = prev;
		prev->cprev = p;

		htable[hash_val] = prev;
		prev->next = entry;

		//if (prev->r1 - prev->cprev->r1 < -128 || prev->r1 - prev->cprev->r1 >127) { cout << "error 4\n"; cout << prev->r1 << "\n" << prev->cprev->r1 << "\n"; }
	}*/


	/*
	* Remove Element at a key/*/

	void pophead(int index)
	{

		HashNode* entry = nullptr;
		entry = htable[index];
		entry = entry->next;
		htable[index] = entry;
		//delete entry;
	}

	//* Search Element at a key

	HashNode* Search(int r, unsigned char character)
	{

		unsigned int hash_val = HashFunc(character);
		HashNode* entry = nullptr;
		entry = htable[hash_val];

		if (entry == NULL)
		{
			entry = nullptr;
			return entry;
		}
		else if (entry->r1 <= r && r <= (entry->r1 + entry->r2 - 1))//was +1 before
		{

			return entry;
		}
		else
		{
			entry = nullptr;
			return entry;
		}
	}

	HashNode* Search(int r, vector<unsigned char> lookahead)
	{
		HashNode* s = nullptr;

		for (std::vector<unsigned char>::iterator it = lookahead.begin(); it != lookahead.end(); ++it)
		{
			s = Search(r, *it);
			if (s != nullptr) { return s; }
		}
		s = nullptr;
		return s;//will be nullptr anyway
	}

	bool modifyr2(HashNode* cposlist)
	{

		if (cposlist->r2 < 127)//127
		{
			cposlist->r2 = cposlist->r2 + 1;
			return false;
		}
		else
		{
			return true; //true means node is full

		}

	}



	int decimalToBinary(int N)
	{

		// To store the binary number 
		ull B_Number = 0;
		int cnt = 0;
		while (N != 0) {
			int rem = N % 2;
			ull c = pow(10, cnt);
			B_Number += rem * c;
			N /= 2;

			// Count used to store exponent value 
			cnt++;
		}

		return B_Number;
	}

	string intToBinaryString(unsigned long n)
	{

		char     bArray[(sizeof(unsigned long) * 8) + 1];
		//index = 32 to fetch each 32 slots
		unsigned index = sizeof(unsigned long) * 8;
		char temp = 0;
		bArray[index] = '\0';
		do {

			temp = (n & 1); // Finding 0 or 1 in LSB
			// Adding ASCII 0 to convert Binary to ASCII
			temp = temp + '0';

			// Copying final ASCII equivalent BIT value (0 or 1)
			bArray[--index] = temp;

		} while (n >>= 1);

		return string(bArray + index);
	}


	string DecimalToBinaryString(int a)
	{
		unsigned int b = (unsigned int)a;
		string binary = "";
		unsigned int mask = 0x80000000u;
		while (mask > 0)
		{
			binary += ((b & mask) == 0) ? '0' : '1';
			mask >>= 1;
		}
		while (binary[0] == '0')
		{
			binary.erase(binary.begin());
		}
		//cout << binary << endl;
		return binary;
	}

	vector<unsigned char> gamma(int val)
	{
		int N = log2(val);
		vector<unsigned char> temp;
		temp.insert(temp.begin(), N, '0');



		string s = DecimalToBinaryString(val);
		//string s = intToBinaryString(val);
		copy(s.begin(), s.end(), back_inserter(temp));

		return temp;
	}



	unsigned int getAbs(int n)
	{
		int const mask = n >> (sizeof(int) * CHAR_BIT - 1);
		return ((n + mask) ^ mask);
	}


	vector<int> degammadepth(vector<unsigned char> gamma)
	{
		int N = 0;
		int x;
		string temp;
		vector<int>output;

		for (int i = 0; i < gamma.size(); i++)
		{
			if (gamma[i] == '0') { N++; }
			else
			{
				x = i + 1;
				temp += gamma[i];
				for (int m = x; m < x + N; m++)
				{
					temp += gamma[m];

				}
				i = x + N - 1;
				//if (N == 0) { i = x; } else{ i = x + N - 1; }
				output.push_back(std::stoi(temp, nullptr, 2) - 1);
				temp.clear();
				N = 0;
				//break; 
			}
		}


		return output;


	}

	vector<int> degammar1(vector<unsigned char> gamma)
	{
		int N = 0;
		int x;
		string temp;
		vector<int>output;

		for (int i = 0; i < gamma.size(); i++)
		{
			if (gamma[i] == '0') { N++; }
			else
			{
				x = i + 1;
				temp += gamma[i];
				for (int m = x; m < x + N; m++)
				{
					temp += gamma[m];

				}
				i = x + N - 1;
				//if (N == 0) { i = x; } else{ i = x + N - 1; }
				if (std::stoi(temp, nullptr, 2) >= 129) { output.push_back(-(std::stoi(temp, nullptr, 2) - 128)); }
				else { output.push_back(std::stoi(temp, nullptr, 2) - 1); }


				temp.clear();
				N = 0;
				//break; 
			}
		}


		return output;


	}

	bool ForwardSearch(HashNode* s, int nextr1, int nextr2)
	{
		s = s->cnext;

		while (s != nullptr)
		{
			if ((s->r1 + s->r2) > (nextr1 + nextr2) && -128 <= (s->r1 - nextr1) && (s->r1 - nextr1) <= 127) //buffer cnode /////////////////////// first inequality was just > before
			{
				return false;
			}
			else if ((s->r1 - nextr1) < -128 || (s->r1 - nextr1) > 127) //other cnode along the way to buffer cnode OR buffer cnode whose r1 delta exceeds 1 bytes capacity
			{
				return true;
			}
			s = s->cnext;

		}
		if (s == nullptr) { return false; }
	}
	vector<unsigned char> output(string::iterator newbitflagbegin) // change to vector<unsigned char> as output
	{
		HashNode* s = head;
		//HashNode *x = nullptr;
		HashNode* p = nullptr;
		//vector<unsigned char> output, buffer;
		vector<unsigned char> r1, c, depth;
		vector<unsigned char>testr2, testr1; // was vector<int>
		int greatestr1 = 0;
		vector<HashNode> hashnode;
		//output.clear();
		//buffer.clear();
		unsigned char newchar;
		unsigned char newr;
		//unsigned char flag = 0;
		//const unsigned char rbit = 1 << 7;
		int counter = 0;
		int sum = 0;
		int i = -1;
		vector<unsigned char> coutput;
		unsigned char newr1;
		int prevr1 = 0;
		int curr1 = 0;
		int redundancy = 0;
		int greatestdepth = 0;
		vector<unsigned char> bin;
		int nextr1 = 0; int nextr2 = 0;
		//int rednextr1 = 0;
		//int rednextr2 = 0;
		//int sub = 0;
		while (s != nullptr)
		{
			/*if (s->r2 == 2)
			{
				*(newbitflagbegin + s->p1) = '0';
				*(newbitflagbegin + s->p2) = '0';
				if (s->r1 + s->r2 - 1 > rednextr1 + rednextr2 - 1)
				{
					sub++;
					rednextr1 = s->r1; //nextr1
					rednextr2 = s->r2; //nextr2
					//sub += s->r1 + s->r2 - 1 - (nextr1 + nextr2 - 1);
				}
				//sub += s->r1 + s->r2 - 1 - (nextr1 + nextr2 - 1);
				s = s->cnext;
				if (s == NULL)
				{

					break;
				}
			}*/

			if (s->r2 == 2 && s->r1 + s->r2 - 1 <= nextr1 + nextr2 - 1) // if rednode is from prevstring
			{
				*(newbitflagbegin + s->p1) = '0';
				*(newbitflagbegin + s->p2) = '0';
				s = s->cnext;
				if (s == NULL)
				{

					break;
				}
			}
			////////////////////////////////////////////////////////////
			/*else if (s->r2 == 2 &&((s->cnext != nullptr&& -128 <= s->cnext->r1 - nextr1 <= 127)|(s->cnext == nullptr))) // rednode is from buffer   // was else if
			{

					*(newbitflagbegin + s->p1) = '0';
					*(newbitflagbegin + s->p2) = '0';
					s = s->cnext;
					if (s == NULL)
					{

						break;
					}


			}*/
			else if (s->r2 == 2 && ForwardSearch(s, nextr1, nextr2) == false)
			{
				*(newbitflagbegin + s->p1) = '0';
				*(newbitflagbegin + s->p2) = '0';
				s = s->cnext;
				if (s == NULL)
				{

					break;
				}
			}
			////////////////////////////////////////////////////////////
			else
			{
				i++;
				if (s->r1 > greatestr1) { greatestr1 = s->r1; }
				if (s->r2 == 2) { redundancy++; } //WAS 0
				else if (s->r2 > greatestdepth) { greatestdepth = s->r2; }
				newchar = s->r2;
				sum += s->r2;
				curr1 = s->r1;

				newr1 = static_cast<unsigned char>(curr1 - prevr1);

				//if (s->r1 > rednextr1 )//+ rednextr2 - 1 //s->r1 > nextr1 + nextr2 - 1
				//{
					//s->r1 -= sub;
				//}


				s->r1 -= nextr1;

				if (s->r1 < -128 || s->r1>127) { cout << "ERROR\n"; cout << "s->r1 delta: " << s->r1 << "\n"; }
				//if (static_cast<int>(static_cast<signed char>(s->r1)) != s->r1) { cout << "CAST ERROR\n"; cout << "s->r1 delta: " << s->r1 << "\n"; }
				hashnode.push_back((HashNode)*s);
				if (s->r1 + s->r2 + nextr1 > nextr1 + nextr2) { nextr1 = s->r1 + nextr1; nextr2 = s->r2; } ///////////////////////////////////////////// was just > before
				//if (s->r1 + s->r2 + nextr1 > rednextr1 + rednextr2) { rednextr1 = s->r1 + nextr1; rednextr2 = s->r2;  } 

				testr2.push_back(s->r2);
				testr1.push_back(s->r1);
				c.push_back(s->c);
				coutput.push_back(s->c);
				coutput.push_back(static_cast<unsigned char>(s->r1)); //was signed
				coutput.push_back(static_cast<unsigned char>(s->r2));


				p = s;
				s = s->cnext;

				if (s == NULL)
				{

					break;
				}
			}



		}
		sum = sum / (i + 1);
		std::cout << "average depth: " << sum << "\n";
		std::cout << "redundant entries: " << redundancy << "\n";
		//std::cout << "greatest depth: " << greatestdepth << " +1 for run length \n";
		std::cout << "greatest depth: " << greatestdepth << "\n";
		std::cout << "greatest r1: " << greatestr1 << "\n";

		cout << "char entries: " << c.size() << "\n";
		cout << "r1 entries: " << testr1.size() << "\n";
		cout << "r2 entries: " << testr2.size() << "\n";
		/*
				string huffc = HuffmanEncode(c);
				string huffr1 = HuffmanEncode(testr1);
				string huffr2 = HuffmanEncode(testr2);
				cout << "Huffman char size: " << huffc.length()/8 << "\n";
				cout << "Huffman r1 size: " << huffr1.length()/8 << "\n";
				cout << "Huffman r2 size: " << huffr2.length()/8 << "\n";
				cout << "Total Huffman coutput size: " << (huffr2.size()/8) + (huffc.size() / 8) + (huffr1.size() / 8) << "\n";
		*/

		cout << "coutput size: " << coutput.size() << "\n";
		//cout << "bin: " << string(bin.begin(),bin.end()) << "\n";
		return coutput;
	}
};



vector<unsigned char> CCdecompress(vector<unsigned char> uoutput, vector<unsigned char> cnodes, vector<unsigned char> bitflag)
{

	//cout << "cnodes size " << cnodes.size() << "\n";
	vector<HashNode> coutput;
	//int ncount = 0;
	for (int i = 0; i < cnodes.size(); i++)
	{
		coutput.push_back(*new HashNode(cnodes[i], static_cast<int>(cnodes[i + 1]), NULL, NULL)); //was unsigned int
		//if (static_cast<int>(cnodes[i + 1]) >128 ) { (coutput.end() - 1)->r1 = 127-static_cast<int>(cnodes[i + 1]); } // WAS 127
		if (static_cast<int>(cnodes[i + 1]) > 127) { (coutput.end() - 1)->r1 = static_cast<int>(cnodes[i + 1]) - 256; } // WAS 127
		//coutput[coutput.size()-1].r2 = static_cast<int>(cnodes[i + 2]);
		/*if ((coutput.end() - 1)->r1 > 125)
		{
			(coutput.end() - 1)->r1 = static_cast<signed int>((coutput.end() - 1)->r1); //newly added in to account for intiial casting of r1 values to unsigned char in output()
		}*/
		(coutput.end() - 1)->r2 = static_cast<int>(cnodes[i + 2]);
		i += 2;
	}
	//cout << "ncount: " << ncount << "\n";

	int maxr = 0;
	int nextr1 = 0;
	int nextr2 = 0;

	for (int l = 0; l < coutput.size(); l++)
	{

		coutput[l].r1 += nextr1;
		//cout << "decomp node \n" << coutput[l].c << " " << coutput[l].r1 << " " << coutput[l].r2 << "\n";
		if (coutput[l].r1 + coutput[l].r2 - 1 > maxr) { maxr = coutput[l].r1 + coutput[l].r2 - 1; }//was 1 since min val of gamma encosing is 1  //WAS+2 on both sides
		if (coutput[l].r1 + coutput[l].r2 > nextr1 + nextr2) //can remove the twos on either side //////////////////////////////////////////////////////////////////////// was just > before
		{
			nextr1 = coutput[l].r1;
			nextr2 = coutput[l].r2;
		}
	}
	//cout << "decomp coutput size " << coutput.size() << "\n";


	vector<vector<unsigned char>> matrix; //reserve size so that random access index can be used. reserved size should be maxr value i.e. greatestr1 + respective node's depth -> can be determined from compressed data
	//maxr = 3;
	//cout << "maxr: " << maxr << "\n";
	//cout << "nextr1+nextr2+2: " << nextr1+nextr2+2 << "\n";
	matrix.resize(maxr);
	vector<unsigned char> decompress;
	vector<unsigned char> chars;
	//int sr1 = 1;
	//int x = 0;
	//vector<HashNode> cinput;
	//undoing r1 delta encoding
	double init = coutput.size();
	progress_bar progress1{ std::clog, 70u, "Decompression Progress: " };
	for (int i = 0; i < coutput.size(); i++)//
	{

		progress1.write((double)i / (double)coutput.size());

		for (int x = 0; x < coutput[i].r2; x++) //was coutput[i].r2 + 2
		{
			matrix[coutput[i].r1 + x - 1].push_back(coutput[i].c);//there was a -1 

		}

	}

	for (int i = 0; i < matrix.size(); i++)
	{

		copy(matrix[i].begin(), matrix[i].end(), back_inserter(chars));
	}

	string newbitflag;
	for (int i = 0; i < bitflag.size(); i++)
	{
		newbitflag += char_to_str(bitflag[i]);
	}

	int c1 = 0; int d1 = 0;
	//cout << "decomp bitflag: " << newbitflag << "\n";
	//cout << "uoutput: " << string(uoutput.begin(), uoutput.end()) << "\n";
	//cout << "chars: " << string(chars.begin(), chars.end()) << "\n";
	for (int i = 0; i < newbitflag.length(); i++)
	{

		if (newbitflag[i] == '0' && d1 < uoutput.size())
		{
			decompress.push_back(uoutput[d1]);

			d1++;


		}
		else if (c1 < chars.size())
		{
			decompress.push_back(chars[c1]);

			c1++;

		}
	}

	return decompress;

}

struct rednode
{
	HashNode* node;
	int p1;
	int p2;
};

// Perform Run Length Encoding (RLE) data compression algorithm
// on string str

void CCencode(vector<unsigned char> input)
{
	std::cout << "CC \n";
	vector<unsigned char> buffer;
	vector<unsigned char> lookahead;
	vector<unsigned char> prevstring;
	//vector<HashNode> coutput;
	vector<unsigned char> coutput;
	bool state;
	//vector<int> positions;
	//positions.clear();
	int counter = 1;
	std::vector<unsigned char>::iterator temp;
	bool atheta = false;
	int i = -1;
	int cposition;
	int newtheta = 0;
	HashNode* cposlist = nullptr;
	HashNode* index = nullptr;
	HashNode* s = nullptr;
	HashNode* prevnode = nullptr;
	HashMap hash;
	int paradox = 0;

	string newbitflag;
	newbitflag.clear();
	newbitflag.append(input.size(), '0');

	int paradoxcount = 0;
	int dollarcount = 0;

	//vector<rednode> rednodes;

	for (std::vector<unsigned char>::iterator it = input.begin(); it != input.end(); ++it)
	{
		++i;

		if ((buffer.empty() == false) && buffer.end() != (temp = find(buffer.begin(), buffer.end(), *it))) //  *it != zerochar && *it != onechar &&       (*it == zerochar || *it == onechar)&&
		{

			cposition = static_cast<int>(distance(buffer.begin(), temp));

			lookahead.clear();

			std::copy(buffer.begin() + cposition + 1, buffer.end(), back_inserter(lookahead));

			prevstring = buffer;
			//cout << string(prevstring.begin(),prevstring.end())<<"\n";
			buffer.clear();
			buffer.push_back(*it);
			++counter;  //was uncommented before
			cposlist = hash.Search(counter - 1, *it);
			//cposlist = hash.Search(counter , *it);

			index = hash.Search(counter - 1, lookahead);

			if (cposlist == nullptr && index == nullptr)
			{
				//rednode node;
				//node.node = 
				hash.Insert(counter - 1, *it, i, i - buffer.size() - (prevstring.size() - cposition) + 1);
				//node.p1 = i;
				//node.p2 = i - buffer.size() - (prevstring.size() - cposition) + 1;
				//rednodes.push_back(node);

				newbitflag[i] = '1';

				newbitflag[i - buffer.size() - (prevstring.size() - cposition) + 1] = '1';

				//if (input[i - buffer.size() - (prevstring.size() - cposition) + 1] != *it) { cout << "ERROR\n"; break; }
			}



			else if (cposlist == nullptr && index != nullptr)
			{
				//rednode node;
				//node.node =	
				hash.InsertBefore(counter - 1, *it, index, i, i - buffer.size() - (prevstring.size() - cposition) + 1);//index wasnt here before
				//node.p1 = i;
				//node.p2 = i - buffer.size() - (prevstring.size() - cposition) + 1;
				//rednodes.push_back(node);

				newbitflag[i] = '1';

				newbitflag[i - buffer.size() - (prevstring.size() - cposition) + 1] = '1';

				//if (input[i - buffer.size() - (prevstring.size() - cposition) + 1] != *it) { cout << "ERROR\n"; break; }
			}

			else
			{

				state = hash.modifyr2(cposlist);

				if (state == false) { newbitflag[i] = '1'; }


			}

			/*for (int x = 0; x < rednodes.size(); x++)
			{
				if (rednodes[x].node->r2 > 2)
				{
					rednodes.erase(rednodes.begin() + x);
					x--;
				}

				else if (rednodes[x].node->r1 + rednodes[x].node->r2 - 1 <= counter - 2) // was else if
				{
					newbitflag[rednodes[x].p1] = '0';
					newbitflag[rednodes[x].p2] = '0';
					hash.pop(rednodes[x].node);
					rednodes.erase(rednodes.begin() + x);
					x--;
				}

			}*/


		}


		else if (prevstring.empty() == false && prevstring.end() != (temp = find(prevstring.begin(), prevstring.end(), *it)))// *it != zerochar && *it != onechar &&         ((*it == zerochar || *it == onechar) &&
		{

			buffer.push_back(*it);

			cposition = static_cast<int>(distance(prevstring.begin(), temp));

			lookahead.clear();
			std::copy(prevstring.begin() + cposition + 1, prevstring.end(), back_inserter(lookahead));


			for (std::vector<unsigned char>::iterator im = lookahead.begin(); im != lookahead.end(); ++im)// to check if common letters occur after cposition and if they also occur before data[i] in the buffer. 
			{
				s = hash.Search(counter - 1, *im);
				if (s != nullptr && find(buffer.begin(), buffer.end(), *im) != buffer.end())
				{
					atheta = true;
					//cout << "trigger 4\n";
					paradox++;
					break;
				}
			}


			index = hash.Search(counter - 1, lookahead);
			cposlist = hash.Search(counter - 1, *it);


			if (atheta == false && (cposlist == nullptr) && index == nullptr)// was s before
			{
				//rednode node;
				//node.node = 
				hash.Insert(counter - 1, *it, i, i - buffer.size() - (prevstring.size() - cposition) + 1);
				//node.p1 = i;
				//node.p2 = i - buffer.size() - (prevstring.size() - cposition) + 1;
				//rednodes.push_back(node);


				newbitflag[i] = '1';

				newbitflag[i - buffer.size() - (prevstring.size() - cposition) + 1] = '1';

				if (input[i - buffer.size() - (prevstring.size() - cposition) + 1] != *it) { cout << "ERROR\n"; break; }
			}

			else if (atheta == false && (cposlist == nullptr) && index != nullptr)// was s or l
			{
				//rednode node;
				//node.node = 
				hash.InsertBefore(counter - 1, *it, index, i, i - buffer.size() - (prevstring.size() - cposition) + 1); //was s or l
				//node.p1 = i;
				//node.p2 = i - buffer.size() - (prevstring.size() - cposition) + 1;
				//rednodes.push_back(node);


				newbitflag[i] = '1';

				newbitflag[i - buffer.size() - (prevstring.size() - cposition) + 1] = '1';

				if (input[i - buffer.size() - (prevstring.size() - cposition) + 1] != *it) { std::cout << "ERROR\n"; break; }
			}

			else if (atheta == false && cposlist != nullptr)
			{
				//cout << "trigger 7\n";
				state = hash.modifyr2(cposlist);

				if (state == false) { newbitflag[i] = '1'; }



			}
			atheta = false;
		}



		else
		{
			//cout << "3rd conditional\n";

			buffer.push_back(*it);

		}


	}

	/*for (int x = 0; x < rednodes.size(); x++) //meant to deal with redundant nodes that exist on counter-1 after compression ends
	{
		if (rednodes[x].node->r2 > 2)
		{
			rednodes.erase(rednodes.begin() + x);
			x--;
		}

		else //else if (rednodes[x].node->r1 + rednodes[x].node->r2 - 1 <= counter - 1) //else if
		{
			newbitflag[rednodes[x].p1] = '0';
			newbitflag[rednodes[x].p2] = '0';
			hash.pop(rednodes[x].node);
			rednodes.erase(rednodes.begin() + x);
			x--;
		}

	}*/

	/*for (int x = 0; x < rednodes.size(); x++)
	{
		if (rednodes[x].node->r2 > 2)
		{
			rednodes.erase(rednodes.begin() + x);
			x--;
		}

		else //if (rednodes[x].node->r1 + rednodes[x].node->r2 - 1 <= counter - 2) // was else if
		{
			newbitflag[rednodes[x].p1] = '0';
			newbitflag[rednodes[x].p2] = '0';
			hash.pop(rednodes[x].node);
			rednodes.erase(rednodes.begin() + x);
			x--;
		}

	}*/
	//cout << "newbitflag " << newbitflag << "\n";
	std::cout << "max r value: " << counter << "\n";
	std::cout << "paradox count: " << paradox << "\n";
	double maxr = counter;
	coutput = hash.output(newbitflag.begin());


	counter = 0;
	vector<unsigned char> uoutput;
	uoutput.clear();


	//	std::ofstream stream("C:/Users/prana/Desktop/uoutput", std::ios::out | std::ios::binary);
	for (int i = 0; i < input.size(); i++)
	{

		if (newbitflag[i] == '0')
		{
			//stream << input[i];
			uoutput.push_back(input[i]);
		}


	}
	//cout << "raw bitfag: " << newbitflag << "\n";
	//std::ifstream in("C:/Users/prana/Desktop/uoutput", std::ios::in | std::ios::binary);
	//std::vector<unsigned char> data((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
	//if (data != uoutput) { cout << "UOUTPUT DOES NOT MATCH\n"; }




	//std::ofstream bitflag("C:/Users/prana/Desktop/bitflag", std::ios::out | std::ios::binary);
	vector<unsigned char> bitflag;
	/*unsigned char bit = 1 << 7;
	int pos=0;
	bitset<8> set; set.reset();
	for (string::iterator i = newbitflag.begin(); i != newbitflag.end(); i++)
	{
		if(*i=='1'){ set[7 - pos] = 1; }

		pos++;
		if (pos == 8)
		{
			pos = 0;
			bitflag.push_back(static_cast<unsigned char>(set.to_ulong()));

		}
	}
	if (newbitflag.length() % 8 != 0) { bitflag.push_back(static_cast<unsigned char>(set.to_ulong())); }*/

	unsigned char flag = 0;
	for (int i = 0; i < newbitflag.size(); i++)
	{
		if (i + 7 < newbitflag.size())//i + 7 <= newbitflag.size()
		{
			flag = static_cast<unsigned char>(stoi(string(newbitflag.begin() + i, newbitflag.begin() + i + 8), nullptr, 2));

			bitflag.push_back(flag);
			i += 7;
			//cout << "trigger\n";
		}
		else
		{
			flag = 0;
			unsigned char bit = 1;
			//cout << "trigger\n";
			//i++;
			//bitset<8> set;

			for (int a = 7; a >= 0; a--)
			{
				if (newbitflag[i] == '1')
				{
					//set[a] = 0;
					flag |= bit << a;
				}
				/*else
				{
					//set[a] = 1;
				}*///cout << "flag: " << char_to_str(flag) << "a: "<<a<<"\n";
				if (i + 1 == newbitflag.length()) { break; }
				i++;
			}
			//bitflag.push_back(static_cast<unsigned char>(set.to_ulong()));
			bitflag.push_back(static_cast<unsigned char>(flag));
			/*string intermediate = string(newbitflag.begin() + i, newbitflag.end());
			cout << "intermediate: " << intermediate << "\n";
			//intermediate.append(8-intermediate.length(),'0');
			//cout << "diff " << 8 - intermediate.length() << "\n";
			//while (intermediate.length() != 8) { intermediate += "0"; }//for (int a = 0; a < (8 - intermediate.length()); a++) { intermediate += "0"; }
			//cout << "intermediate: " << intermediate << "\n";
			flag = static_cast<unsigned char>(stoi(intermediate), nullptr, 2);

			flag = flag<< 8 - intermediate.length();
			cout << "bin dtr of flag: " << char_to_str(flag) << "\n";
			cout << "int " << static_cast<int>(flag) << "\n";*/
			//bitflag.push_back(flag);
			break;
		}

	}


	//std::ifstream ii("C:/Users/prana/Desktop/bitflag", std::ios::in | std::ios::binary);
	//std::vector<unsigned char> b((std::istreambuf_iterator<char>(ii)), std::istreambuf_iterator<char>());
	//string inp; inp.assign((std::istreambuf_iterator<char>(ii)), (std::istreambuf_iterator<char>()));
	//cout << "b size in bytes: " << b.size() << "\n"; // inp.length()
	//cout << "newbitflag bits: " << newbitflag.length() << "\n";




	//std::cout << "CC of coutput: " << coutput.size() << "\n";
	std::cout << "CC of uoutput: " << uoutput.size() << "\n";

	//std::cout << "CC of newbitflag: " << newbitflag.length() << " bits\n";
	cout << "bitflag bytes: " << bitflag.size() << "\n";

	//PackBits(newbitflag);

	cout << "total compressed size: " << bitflag.size() + coutput.size() + uoutput.size() << "\n";

	//recursive CC
	/*if (bitflag.size() + coutput.size() + uoutput.size() < input.size())
	{
		CCencode(uoutput);
	}
	else { cout << "STOP\n"; }*/

	//HUFFMAN OF UOUTPUT
	/*int ufrequencies[UniqueSymbols] = { 0 };
	string ubinaryoutput;
	if (uoutput.empty() == false)
	{

		for (std::vector<unsigned char>::iterator it = uoutput.begin(); it != uoutput.end(); ++it)
		{
			++ufrequencies[static_cast<int>(*it)];
		}
		INode* uroot = BuildTree(ufrequencies);

		HuffCodeMap ucodes;
		GenerateCodes(uroot, HuffCode(), ucodes);
		delete uroot;

		for (std::vector<unsigned char>::iterator it = uoutput.begin(); it != uoutput.end(); ++it)
		{
			ubinaryoutput += ucodes.find(*it)->second;
		}
	}
	cout << "Huffman uoutput size: " << ubinaryoutput.length() / 8 << "\n";


	vector<unsigned char> delta;
	int previndex = 0;
	for (int i = 0; i < newbitflag.length(); i++)
	{
		if (newbitflag[i] == '1')
		{
			delta.push_back(static_cast<unsigned char>(i - previndex));
			previndex = i;
		}
	}
	string deltapos = HuffmanEncode(delta);
	cout << "delta size: " << delta.size()<< "\n";
	cout << "delta huffman size: " << deltapos.length() / 8 << "\n";
	*/
	vector<unsigned char> decompress = CCdecompress(uoutput, coutput, bitflag);

	cout << "\ncomparison in progress: \n";
	if (decompress == input)
	{
		cout << "DECOMPRESSION SUCCESS\n";
	}
	else {
		cout << "FAIL\n" << "\n"; //cout<<string(decompress.begin(), decompress.end()) << "\n";
	}

	//return uoutput, coutput, bitflag;
	//return uoutput;






}

/*
vector<int> degammadepth(vector<unsigned char> gamma)
{
	int N = 0;
	int x;
	string temp;
	vector<int>output;

	for (int i = 0; i < gamma.size(); i++)
	{
		if (gamma[i] == '0') { N++; }
		else
		{
			x = i + 1;
			temp += gamma[i];
			for (int m = x; m < x + N; m++)
			{
				temp += gamma[m];

			}
			i = x + N - 1;
			//if (N == 0) { i = x; } else{ i = x + N - 1; }
			output.push_back(std::stoi(temp, nullptr, 2) - 1);
			temp.clear();
			N = 0;
			//break;
		}
	}


	return output;


}

vector<int> degammar1(vector<unsigned char> gamma)
{
	int N = 0;
	int x;
	string temp;
	vector<int>output;

	for (int i = 0; i < gamma.size(); i++)
	{
		if (gamma[i] == '0') { N++; }
		else
		{
			x = i + 1;
			temp += gamma[i];
			for (int m = x; m < x + N; m++)
			{
				temp += gamma[m];

			}
			i = x + N - 1;
			//if (N == 0) { i = x; } else{ i = x + N - 1; }
			if (std::stoi(temp, nullptr, 2) >= 129) { output.push_back(-(std::stoi(temp, nullptr, 2) - 128)); }
			else { output.push_back(std::stoi(temp, nullptr, 2) - 1); }


			temp.clear();
			N = 0;
			//break;
		}
	}


	return output;


}*/

string PackBitsDecompress(vector<unsigned char> input)
{

	string bin;
	string output;
	for (int i = 0; i < input.size(); i++)
	{
		bin += char_to_str(input[i]);
	}

	for (string::iterator i = bin.begin(); i != bin.end(); i++)
	{
		if (distance(bin.begin(), i) + 7 > bin.length() - 1) { break; }
		//cout << "distance: " << distance(bin.begin(), i)<<"\n";
		if (*i == '1')
		{
			++i;
			int run = stoi(string(i, i + 8), nullptr, 2) & 255;
			//run++;
			run += 2;
			//cout << "run: " << run << "\n";
			i += 8;
			unsigned char c = stoi(string(i, i + 8), nullptr, 2) & 255;
			i += 7;

			output.append(run, c);

		}
		else
		{
			++i;
			//if (distance(bin.begin(),i) + 7 > bin.length() - 1) { break; }
			int run = stoi(string(i, i + 8), nullptr, 2) & 255;
			//cout << "run: " << run << "\n";
			//run++;
			run += 1;///////////
			i += 8;
			for (int x = run; x > 0; x--)
			{
				if (distance(bin.begin(), i) + 7 > bin.length() - 1) { break; }
				unsigned char c = stoi(string(i, i + 8), nullptr, 2) & 255;
				i += 8;
				output += c;
			}
			if (distance(bin.begin(), i) + 7 > bin.length() - 1) { break; }
			i--;
		}
	}
	//cout << "output: " << output << "\n";
	return output;

}



void PackBits(string input)
{
	// Open files

	string out;

	int c, c1 = -1;  // current and last char
	char t[256];  // buffer for literals
	int run = 0;  // current run length
	int lrun = 0;  // current literal run length
	int i = -1;
	// Compress

		// RLE encode
	do {
		++i;
		//cout << "i: " << i << "\n";
		c = static_cast<int>(input[i]);
		//cout << "i: " << i << "\n";
		if (c == c1) {
			if (lrun > 1)
			{
				out += '0';
				out += char_to_str(static_cast<unsigned char>(lrun - 2));
				//cout << "lrun: " << lrun << "\n";
				for (int x = 0; x < lrun - 1; x++)///////////////
				{
					out += char_to_str(static_cast<unsigned char>(t[x]));
					//cout << "t[x]: " << t[x] << "\n";
				}
				//fwrite(t, 1, lrun - 1, out.end()); 
			}
			lrun = 0;
			if (++run == 256) { out += '1'; out += char_to_str(static_cast<unsigned char>(255)); out += char_to_str(static_cast<unsigned char> (c1)); run = 0; t[lrun++] = c; }
		}
		else
		{
			if (run > 0) { out += '1'; out += char_to_str(static_cast<unsigned char>(run - 1)); out += char_to_str(static_cast<unsigned char> (c1)); run = 0; }
			if (lrun == 256)
			{
				out += '0';
				out += char_to_str(static_cast<unsigned char>(255));
				for (int x = 0; x < 256; x++)
				{
					out += char_to_str(static_cast<unsigned char>(t[x]));
				}
				//fwrite(t, 1, 256, out); 
				lrun = 0;
			}
			t[lrun++] = c;
		}
		c1 = c;
	} while (i <= input.length() - 1); //c != EOF
	if (lrun > 1)
	{
		out += '0';
		out += char_to_str(static_cast<unsigned char>(lrun - 2));
		for (int x = 0; x < lrun - 1; x++)
		{
			out += char_to_str(static_cast<unsigned char>(t[x]));
			//cout << "t[x]: " << t[x] << "\n";
		}
		//fwrite(t, 1, lrun - 1, out); 

	}
	vector<unsigned char>output;
	/*for (string::iterator it = out.begin(); it != out.end(); it++)
	{
		if (distance(out.begin(),it)+7<=out.size()-1)
		{
			int l = stoi(string(it, it + 8), nullptr, 2);
			output.push_back(static_cast<unsigned char>(l));
			it += 6;
		}
		else
		{
			int l = stoi(string(it, out.end()), nullptr, 2);
			output.push_back(static_cast<unsigned char>(l));
			break;
		}
	}*/

	unsigned char flag = 0;
	for (int i = 0; i < out.size(); i++)
	{
		if (i + 7 < out.size())//i + 7 <= newbitflag.size()
		{
			flag = static_cast<unsigned char>(stoi(string(out.begin() + i, out.begin() + i + 8), nullptr, 2));

			output.push_back(flag);
			i += 7;
			//cout << "trigger\n";
		}
		else
		{
			flag = 0;
			unsigned char bit = 1;


			for (int a = 7; a >= 0; a--)
			{
				if (out[i] == '1')
				{

					flag |= bit << a;
				}

				if (i + 1 == out.length()) { break; }
				i++;
			}

			output.push_back(static_cast<unsigned char>(flag));

			break;
		}

	}
	cout << "Packbits of 9 bit byte size: " << output.size() << "\n";
	//cout << "Packbits size: " << out.length()/8 << "\n";
	//cout << "out: " << out << "\n";
	string decompress = PackBitsDecompress(output);
	if (decompress == input) { cout << "PackBits decompression Success\n"; }
	else { cout << "PackBits decompression FAIL\n"; }
}

int main()
{
	//string data;
	//std::ifstream ifs("C:/Users/prana/Desktop/enwik8", std::ios::binary);
	//example_static();
	//std::ifstream ifs("E:/pranavWork/comrpession project/compression/silesia/dickens", std::ios::binary);
	//std::ifstream stream("C:/Users/prana/Desktop/Pranav_Venkatram_passport.jpg", std::ios::binary);
	//std::ifstream ifs("C:/Users/prana/Downloads/rgb8bit/artificial.ppm", std::ios::binary);
	//std::ifstream stream("C:/Users/prana/Downloads/rgb8bit/flower_foveon.ppm", std::ios::binary);
	//std::ifstream stream("E:/pranavWork/comrpession project/DNA/chmpxx", std::ios::binary);

	//data.assign((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));

	std::ifstream stream("D:/CC/silesia/dickens", std::ios::in | std::ios::binary);
	//std::ifstream stream("E:/pranavWork/comrpession_project/compression/pizza_and_chili_repetitive_corpus/Real/world_leaders", std::ios::in | std::ios::binary);
	//std::ifstream stream("E:/pranavWork/comrpession_project/compression/1034.db", std::ios::in | std::ios::binary);

	std::vector<unsigned char> data((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());

	//string test = "THEPHONEBLAH";
	//string test = "HNAEHSEHSES";
	//vector<unsigned char> data(test.begin(), test.end());
	//std::ifstream stream("C:/Users/prana/Desktop/enwik8", std::ios::in | std::ios::binary);
	//vector<unsigned char>data(i.begin()+3731830, i.begin() + 3732290);
	cout << "open file\n";
	//string t = "ABABACXDBCDX";
	//string t = "ALMBNOAPQBRSATUCVWXYDBZCDX";
	//string t = "This file contains the 'main' function. Program execution begins and ends there.";
	//string t = "THEPHONEBLAH";
	//vector<unsigned char> data(t.begin(), t.end());
	cout << "input length: " << data.size() << "\n";
	//cout<<"assign to string\n";

	//lz77::compress_t lzcompress(8, 4096);
	//std::string lzcompressed = lzcompress.feed(data);
	//cout << "LZ77 size:  " << lzcompressed.length() << "\n";


	//PackBits(string(data.begin(), data.end()));//////////////////////////////////////////////////////////////////

	//string l; l.append(256, 'A');
	//PackBits(l);

	CCencode(data);

	/*
	vector<unsigned char> decompress = CCdecompress(uoutput,coutput,bitflag);

	cout << "\ncomparison in progress: \n";
	if (decompress == data)
	{
		cout << "DECOMPRESSION SUCCESS\n";
	}
	else {
		cout << "FAIL\n";
	}*/


	return 0;
}