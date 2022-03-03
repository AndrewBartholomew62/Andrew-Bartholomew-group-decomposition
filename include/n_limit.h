class n_limit;

typedef pair<int,n_limit*> threshold;

class n_limit
{

public:

	int n;
	list<threshold> thresholds;

	n_limit(vector<int> z);
	n_limit(const n_limit& l);
	~n_limit();
	
	
	void print(ostream& oss, int indent, string prefix);

};

ostream& operator << (ostream& os, threshold t);
ostream& operator << (ostream& os, list<threshold> thresholds);
bool operator > (n_limit& l, vector<int> u);
bool operator == (n_limit&l1, n_limit& l2);
bool operator != (n_limit&l1, n_limit& l2);
