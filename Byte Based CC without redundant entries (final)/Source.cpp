#include <iostream>
#include <vector>
#include <fstream>
#include <bitset>
#include <string>
#include <iomanip>

const int TABLE_SIZE = 256;

using namespace std;

class progress_bar
{
	static const auto overhead = sizeof " [100%]";

	std::ostream &os;
	const std::size_t bar_width;
	std::string message;
	const std::string full_bar;

public:
	progress_bar(std::ostream &os, std::size_t line_width,
				 std::string message_, const char symbol = '.')
		: os{os},
		  bar_width{line_width - overhead},
		  message{std::move(message_)},
		  full_bar{std::string(bar_width, symbol) + std::string(bar_width, ' ')}
	{
		if (message.size() + 1 >= bar_width || message.find('\n') != message.npos)
		{
			os << message << '\n';
			message.clear();
		}
		else
		{
			message += ' ';
		}
		write(0.0);
	}

	// not copyable
	progress_bar(const progress_bar &) = delete;
	progress_bar &operator=(const progress_bar &) = delete;

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

string char_to_str(unsigned char ch)
{
	bitset<8> b(ch);
	return b.to_string();
}

// common letter node in hash table declaration
class HashNode
{
public:
	unsigned char c; //key is same as char
	int r1, r2, p1, p2;
	HashNode *next;
	HashNode *cnext;
	HashNode *cprev;
	HashNode(unsigned char c, int r1, int p1, int p2)
	{
		this->c = c;
		this->r1 = r1;
		this->r2 = 2;
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
	HashNode **htable;
	HashNode *head, *tail;

public:
	HashMap()
	{
		head = NULL;
		tail = NULL;
		htable = new HashNode *[TABLE_SIZE];
		for (int i = 0; i < TABLE_SIZE; i++)
			htable[i] = NULL;
	}
	~HashMap()
	{
		for (int i = 0; i < TABLE_SIZE; ++i)
		{
			HashNode *entry = htable[i];
			while (entry != NULL)
			{
				HashNode *prev = entry;
				entry = entry->next;
				delete prev;
			}
		}
		delete[] htable;
	}
	/*
	* Hash Function
	*/
	int HashFunc(unsigned char c)
	{
		int key = static_cast<int>(c);
		return key & 255;
	}

	/*
	* Insert Element at a key at the start of the list
	*/
	void Insert(int r1, unsigned char c, int p1, int p2)
	{
		unsigned int hash_val = HashFunc(c);
		HashNode *prev = nullptr;
		HashNode *p = nullptr;
		HashNode *entry = nullptr;
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
	}

	void InsertBefore(int r1, unsigned char c, HashNode *index, int p1, int p2)
	{

		unsigned int hash_val = HashFunc(c);
		HashNode *prev = nullptr;
		HashNode *p = nullptr;
		HashNode *entry = nullptr;
		entry = htable[hash_val];
		prev = new HashNode(c, r1, p1, p2);

		if (index->cprev != nullptr)
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
	}

	/*
	* Remove Element at a key/*/

	void pophead(int index)
	{

		HashNode *entry = nullptr;
		entry = htable[index];
		entry = entry->next;
		htable[index] = entry;
	}

	//* Search Element at a key

	HashNode *Search(int r, unsigned char character)
	{

		unsigned int hash_val = HashFunc(character);
		HashNode *entry = nullptr;
		entry = htable[hash_val];

		if (entry == NULL)
		{
			entry = nullptr;
			return entry;
		}
		else if (entry->r1 <= r && r <= (entry->r1 + entry->r2 - 1))
		{

			return entry;
		}
		else
		{
			entry = nullptr;
			return entry;
		}
	}

	HashNode *Search(int r, vector<unsigned char> lookahead)
	{
		HashNode *s = nullptr;

		for (std::vector<unsigned char>::iterator it = lookahead.begin(); it != lookahead.end(); ++it)
		{
			s = Search(r, *it);
			if (s != nullptr)
			{
				return s;
			}
		}
		s = nullptr;
		return s;
	}

	bool modifyr2(HashNode *cposlist)
	{

		if (cposlist->r2 < 127)
		{
			cposlist->r2 = cposlist->r2 + 1;
			return false;
		}
		else
		{
			return true;
		}
	}

	bool ForwardSearch(HashNode *s, int nextr1, int nextr2)
	{
		s = s->cnext;

		while (s != nullptr)
		{
			if ((s->r1 + s->r2) > (nextr1 + nextr2) && -128 <= (s->r1 - nextr1) && (s->r1 - nextr1) <= 127)
			{
				return false;
			}
			else if ((s->r1 - nextr1) < -128 || (s->r1 - nextr1) > 127)
			{
				return true;
			}
			s = s->cnext;
		}
		if (s == nullptr)
		{
			return false;
		}
	}

	vector<unsigned char> output(string::iterator newbitflagbegin)
	{
		HashNode *s = head;
		HashNode *p = nullptr;
		vector<unsigned char> r1, c, depth;
		vector<unsigned char> testr2, testr1;
		int greatestr1 = 0;
		vector<HashNode> hashnode;
		unsigned char newchar;
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
		int nextr1 = 0;
		int nextr2 = 0;

		while (s != nullptr)
		{
			if (s->r2 == 2 && s->r1 + s->r2 - 1 <= nextr1 + nextr2 - 1)
			{
				*(newbitflagbegin + s->p1) = '0';
				*(newbitflagbegin + s->p2) = '0';
				s = s->cnext;
				if (s == NULL)
				{

					break;
				}
			}

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
			else
			{
				i++;
				if (s->r1 > greatestr1)
				{
					greatestr1 = s->r1;
				}
				if (s->r2 == 2)
				{
					redundancy++;
				}
				else if (s->r2 > greatestdepth)
				{
					greatestdepth = s->r2;
				}
				newchar = s->r2;
				sum += s->r2;
				curr1 = s->r1;

				newr1 = static_cast<unsigned char>(curr1 - prevr1);
				s->r1 -= nextr1;

				if (s->r1 < -128 || s->r1 > 127)
				{
					cout << "ERROR\n";
					cout << "s->r1 delta: " << s->r1 << "\n";
				}
				hashnode.push_back((HashNode)*s);
				if (s->r1 + s->r2 + nextr1 > nextr1 + nextr2)
				{
					nextr1 = s->r1 + nextr1;
					nextr2 = s->r2;
				}

				testr2.push_back(s->r2);
				testr1.push_back(s->r1);
				c.push_back(s->c);
				coutput.push_back(s->c);
				coutput.push_back(static_cast<unsigned char>(s->r1));
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
		std::cout << "greatest depth: " << greatestdepth << "\n";
		std::cout << "greatest r1: " << greatestr1 << "\n";

		cout << "char entries: " << c.size() << "\n";
		cout << "r1 entries: " << testr1.size() << "\n";
		cout << "r2 entries: " << testr2.size() << "\n";

		cout << "coutput size: " << coutput.size() << "\n";
		return coutput;
	}
};

vector<unsigned char> CCdecompress(vector<unsigned char> uoutput, vector<unsigned char> cnodes, vector<unsigned char> bitflag)
{
	vector<HashNode> coutput;
	for (int i = 0; i < cnodes.size(); i++)
	{
		coutput.push_back(*new HashNode(cnodes[i], static_cast<int>(cnodes[i + 1]), NULL, NULL));

		if (static_cast<int>(cnodes[i + 1]) > 127)
		{
			(coutput.end() - 1)->r1 = static_cast<int>(cnodes[i + 1]) - 256;
		}

		(coutput.end() - 1)->r2 = static_cast<int>(cnodes[i + 2]);
		i += 2;
	}

	int maxr = 0;
	int nextr1 = 0;
	int nextr2 = 0;

	for (int l = 0; l < coutput.size(); l++)
	{

		coutput[l].r1 += nextr1;
		if (coutput[l].r1 + coutput[l].r2 - 1 > maxr)
		{
			maxr = coutput[l].r1 + coutput[l].r2 - 1;
		}
		if (coutput[l].r1 + coutput[l].r2 > nextr1 + nextr2)
		{
			nextr1 = coutput[l].r1;
			nextr2 = coutput[l].r2;
		}
	}

	vector<vector<unsigned char>> matrix;
	matrix.resize(maxr);
	vector<unsigned char> decompress;
	vector<unsigned char> chars;

	double init = (double)coutput.size();
	progress_bar progress1{std::clog, 70u, "Decompression Progress: "};
	for (int i = 0; i < coutput.size(); i++)
	{

		progress1.write((double)i / (double)coutput.size());

		for (int x = 0; x < coutput[i].r2; x++)
		{
			matrix[coutput[i].r1 + x - 1].push_back(coutput[i].c);
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

	int c1 = 0;
	int d1 = 0;

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
	HashNode *node;
	int p1;
	int p2;
};

void CCencode(vector<unsigned char> input)
{
	std::cout << "CC \n";
	vector<unsigned char> buffer;
	vector<unsigned char> lookahead;
	vector<unsigned char> prevstring;
	vector<unsigned char> coutput;
	bool state;
	int counter = 1;
	std::vector<unsigned char>::iterator temp;
	bool atheta = false;
	int i = -1;
	int cposition;
	int newtheta = 0;
	HashNode *cposlist = nullptr;
	HashNode *index = nullptr;
	HashNode *s = nullptr;
	HashNode *prevnode = nullptr;
	HashMap hash;
	int paradox = 0;

	string newbitflag;
	newbitflag.clear();
	newbitflag.append(input.size(), '0');

	int paradoxcount = 0;
	int dollarcount = 0;

	for (std::vector<unsigned char>::iterator it = input.begin(); it != input.end(); ++it)
	{
		++i;

		if ((buffer.empty() == false) && buffer.end() != (temp = find(buffer.begin(), buffer.end(), *it)))
		{

			cposition = static_cast<int>(distance(buffer.begin(), temp));

			lookahead.clear();

			std::copy(buffer.begin() + cposition + 1, buffer.end(), back_inserter(lookahead));

			prevstring = buffer;
			buffer.clear();
			buffer.push_back(*it);
			++counter;
			cposlist = hash.Search(counter - 1, *it);

			index = hash.Search(counter - 1, lookahead);

			if (cposlist == nullptr && index == nullptr)
			{
				hash.Insert(counter - 1, *it, i, i - buffer.size() - (prevstring.size() - cposition) + 1);

				newbitflag[i] = '1';

				newbitflag[i - buffer.size() - (prevstring.size() - cposition) + 1] = '1';
			}

			else if (cposlist == nullptr && index != nullptr)
			{
				hash.InsertBefore(counter - 1, *it, index, i, i - buffer.size() - (prevstring.size() - cposition) + 1);

				newbitflag[i] = '1';

				newbitflag[i - buffer.size() - (prevstring.size() - cposition) + 1] = '1';
			}

			else
			{

				state = hash.modifyr2(cposlist);

				if (state == false)
				{
					newbitflag[i] = '1';
				}
			}
		}

		else if (prevstring.empty() == false && prevstring.end() != (temp = find(prevstring.begin(), prevstring.end(), *it)))
		{

			buffer.push_back(*it);

			cposition = static_cast<int>(distance(prevstring.begin(), temp));

			lookahead.clear();
			std::copy(prevstring.begin() + cposition + 1, prevstring.end(), back_inserter(lookahead));

			for (std::vector<unsigned char>::iterator im = lookahead.begin(); im != lookahead.end(); ++im)
			{
				s = hash.Search(counter - 1, *im);
				if (s != nullptr && find(buffer.begin(), buffer.end(), *im) != buffer.end())
				{
					atheta = true;
					paradox++;
					break;
				}
			}

			index = hash.Search(counter - 1, lookahead);
			cposlist = hash.Search(counter - 1, *it);

			if (atheta == false && (cposlist == nullptr) && index == nullptr)
			{

				hash.Insert(counter - 1, *it, i, i - buffer.size() - (prevstring.size() - cposition) + 1);

				newbitflag[i] = '1';

				newbitflag[i - buffer.size() - (prevstring.size() - cposition) + 1] = '1';

				if (input[i - buffer.size() - (prevstring.size() - cposition) + 1] != *it)
				{
					cout << "ERROR\n";
					break;
				}
			}

			else if (atheta == false && (cposlist == nullptr) && index != nullptr)
			{

				hash.InsertBefore(counter - 1, *it, index, i, i - buffer.size() - (prevstring.size() - cposition) + 1);

				newbitflag[i] = '1';

				newbitflag[i - buffer.size() - (prevstring.size() - cposition) + 1] = '1';

				if (input[i - buffer.size() - (prevstring.size() - cposition) + 1] != *it)
				{
					std::cout << "ERROR\n";
					break;
				}
			}

			else if (atheta == false && cposlist != nullptr)
			{
				state = hash.modifyr2(cposlist);

				if (state == false)
				{
					newbitflag[i] = '1';
				}
			}
			atheta = false;
		}

		else
		{

			buffer.push_back(*it);
		}
	}

	std::cout << "max r value: " << counter << "\n";
	std::cout << "paradox count: " << paradox << "\n";
	double maxr = counter;
	coutput = hash.output(newbitflag.begin());

	counter = 0;
	vector<unsigned char> uoutput;
	uoutput.clear();

	for (int i = 0; i < input.size(); i++)
	{

		if (newbitflag[i] == '0')
		{
			uoutput.push_back(input[i]);
		}
	}

	vector<unsigned char> bitflag;

	unsigned char flag = 0;
	for (int i = 0; i < newbitflag.size(); i++)
	{
		if (i + 7 < newbitflag.size())
		{
			flag = static_cast<unsigned char>(stoi(string(newbitflag.begin() + i, newbitflag.begin() + i + 8), nullptr, 2));

			bitflag.push_back(flag);
			i += 7;
		}
		else
		{
			flag = 0;
			unsigned char bit = 1;

			for (int a = 7; a >= 0; a--)
			{
				if (newbitflag[i] == '1')
				{
					flag |= bit << a;
				}

				if (i + 1 == newbitflag.length())
				{
					break;
				}
				i++;
			}

			bitflag.push_back(static_cast<unsigned char>(flag));

			break;
		}
	}

	std::cout << "CC of uoutput: " << uoutput.size() << "\n";

	cout << "bitflag bytes: " << bitflag.size() << "\n";

	cout << "total compressed size: " << bitflag.size() + coutput.size() + uoutput.size() << "\n";

	vector<unsigned char> decompress = CCdecompress(uoutput, coutput, bitflag);

	cout << "\ncomparison in progress: \n";
	if (decompress == input)
	{
		cout << "DECOMPRESSION SUCCESS\n";
	}
	else
	{
		cout << "FAIL\n"
			 << "\n";
	}
}

int main()
{

	std::ifstream stream("D:/Everything else/CC/silesia/dickens", std::ios::in | std::ios::binary);

	std::vector<unsigned char> data((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());

	cout << "open file\n";

	cout << "input length: " << data.size() << "\n";

	CCencode(data);

	return 0;
}