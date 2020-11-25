#define NOMINMAX
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <queue>
#include <stack>
#include <Windows.h>

using namespace std;

#pragma region defines
#define F first
#define S second
#define DOWN 0
#define RIGHT 1
#define UP 2
#define LEFT 3
#define debug(x) cout << #x << ": " << x << '\n'
#define BLOCK '!'
#define WATER '.'
typedef long long ll;
typedef pair<int, int> pii;
#pragma endregion

#pragma region global var
int ar[] = { 1, 0, -1, 0 };
int ac[] = { 0, 1, 0, -1 };

//constants at any time in the game
int n, m, k;
vector<vector<class cell>> cells;
vector<vector<int>> bounds[4];

//change with game state
char grid[1005][1005];
int color[1005][1005];
bool vis[1005 * 1005];
#pragma endregion

#pragma region Building data structuer
#pragma region CELL
namespace std {
	template <>
	struct hash<pii> {
		size_t operator()(const pii& p) const {
			hash<int> hasher;
			return ((hasher(p.first) ^ (hasher(p.second) << 1)) >> 1);
		}
	};
}

class cell : public pii {
public:

	static bool valid(cell c, int n, int m) {
		return c.first < n && c.second < m && c.first >= 0 && c.second >= 0;
	}
	
	cell() : pii() {};

	cell(const int& x, const int& y) : pii(x, y) {}

	cell(const pii& p) : pii(p) {}

	cell getDown() {
		return cell(first + ar[DOWN], second + ac[DOWN]);
	}

	cell getRight() {
		return cell(first + ar[RIGHT], second + ac[RIGHT]);

	}

	cell getUp() {
		return cell(first + ar[UP], second + ac[UP]);
	}

	cell getLeft() {
		return cell(first + ar[LEFT], second + ac[LEFT]);
	}

	cell dirNeighbour(int dir) {
		switch (dir) {
		case DOWN:
			return getDown();
		case RIGHT:
			return getRight();
		case UP:
			return getUp();
		case LEFT:
			return getLeft();
		default:
			break;
		}
		return cell(-1, -1);
	}

	vector<cell> neighbours() {
		vector<cell> ve;
		ve.resize(4);
		ve[DOWN] = getDown();
		ve[UP] = getUp();
		ve[RIGHT] = getRight();
		ve[LEFT] = getLeft();
		return ve;
	}

	cell operator +(const pii& other) const {
		return cell(first + other.first, second + other.second);
	}

	cell operator -(const pii& other) const {
		return cell(first - other.first, second - other.second);
	}

	friend cell abs(const cell& c);

	static int manhaten(const cell& c1, const cell& c2) {
		cell res = abs(c1 - c2);
		return res.first + res.second;
	}

	friend ostream &operator<<(ostream &os, const cell  &c);

	friend struct hash<cell>;
};

cell abs(const cell& c) {
	return cell(abs(c.first), abs(c.second));
}

namespace std {
	template <>
	struct hash<cell> {
		size_t operator()(const cell& cell) const {
			hash<int> hasher;
			return ((hasher(cell.first) ^ (hasher(cell.second) << 1)) >> 1);
			//const int prime = 31;
			//int res = 1;
			//res = res * prime + cell.first;
			//res = res * prime + cell.second;
			//return res;
		}
	};
}

ostream &operator<<(ostream &os, const cell& c) {
	return os << '(' << c.first << ", " << c.second << ')';
}

bool ok(const cell& x, const cell& y) {
	if (x.second == y.second)
		if (x.first < y.first)
			return true;
		else
			return false;
	else
		if (x.second < y.second)
			return true;
		else
			return false;
}
#pragma endregion

#pragma region PIECE
class piece {
public:
	pii shifts;
	int color;
	piece() : color(), shifts() {};
	piece(const int& color, const pii shifts = pii(0, 0)) : shifts(shifts), color(color) {};

	int validMove(int dir, int amount = 1) const {
		set<int> indices;
		for (auto idx : bounds[dir][color])
			indices.insert(idx);

		for (int i = 1; i <= amount; i++)
			for (auto it = indices.begin(); it != indices.end(); it++) {
				cell nc = cells[color][*it] + shifts + pii(i * ar[dir], i * ac[dir]);
				if (cell::valid(nc, n, m))
					if (grid[nc.first][nc.second] == WATER)
						continue;
					else
						if (::color[nc.first][nc.second] == color) {
							auto it1 = it;
							++it;
							indices.erase(it1);
							--it;
							continue;
						}
						else
							return i - 1;
				else
					return i - 1;
			}
		return amount;
	}

	size_t size() {
		return cells[color].size();
	}

	virtual cell operator [] (int idx) const {
		return cells[color][idx] + shifts;
	}

	virtual bool operator ==(const piece& other) const {
		return shifts == other.shifts && color == other.color && getChar() == other.getChar();
	}

	virtual bool operator !=(const piece& other) const {
		return !(*this == other);
	}

	virtual char getChar() const = 0;

	virtual vector<piece*> nextMoves() const = 0;

	friend ostream &operator<<(ostream &os, const piece &p);

	friend struct hash<piece*>;
};

namespace std {
	template <>
	struct hash<piece*> {
		size_t operator()(const piece* p) const {
			hash<int> intHasher;
			hash<pii> pairHasher;
			return ((intHasher(p->color) ^ (pairHasher(p->shifts) << 1)) >> 1);
			const int prime = 31;
			int res = 1;
			res = res * prime + p->shifts.first;
			res = res * prime + p->shifts.second;
			res = res * prime + p->color;
			return res;
		}
	};
}

ostream &operator<<(ostream &os, const piece& p) {
	os << "[";
	for (const auto& c : cells[p.color])
		os << c + p.shifts;
	os << "]";
	return os;
}

class bridge : public piece {
public:
	static const char CHAR = '*';
	bridge(const int& color) : piece(color) {};
	bridge(const int& color, const pii& shifts) : piece(color, shifts) {};

	virtual char getChar() const {
		return CHAR;
	}

	virtual vector<piece*> nextMoves() const {
		vector<piece*> res;
		for (int dir = 0; dir < 4; dir++) {
			//int amount = validMove(dir, max(n, m));
			//for (int i = 1; i <= amount; i++)
			//	res.push_back(new bridge(color, pii(shifts.first + ar[dir] * i, shifts.second + ac[dir] * i)));
			if (validMove(dir))
				res.push_back(new bridge(color, pii(shifts.first + ar[dir], shifts.second + ac[dir])));
		}
		return res;
	}

};

class grass : public piece {
public:
	static const char CHAR = '#';
	grass(const int& color) : piece(color) {};
	grass(const int& color, const pii& shifts) : piece(color, shifts) {};

	virtual char getChar() const {
		return CHAR;
	}

	virtual vector<piece*> nextMoves() const {
		vector<piece*> res;
		for (int dir = 0; dir < 4; dir++) {
			//int amount = validMove(dir, max(n, m));
			//for (int i = 1; i <= amount; i++)
			//	res.push_back(new grass(color, pii(shifts.first + ar[dir] * i, shifts.second + ac[dir] * i)));
			if (validMove(dir))
				res.push_back(new grass(color, pii(shifts.first + ar[dir], shifts.second + ac[dir])));
		}

		return res;
	}

};
#pragma endregion

#pragma region STATE
class state {
public:
	state() : pieces() {}

	state(const vector<piece*>& pieces) : pieces(pieces) {}

	/*state(const state& s) {
		for (const auto& p : s.pieces)
			if (p->getChar() == Wpiece::CHAR)
				pieces.push_back(new Wpiece(*dynamic_cast<Wpiece*>(p)));
			else
				if(p->getChar() == Gpiece::CHAR)
					pieces.push_back(new Gpiece(*dynamic_cast<Gpiece*>(p)));
	}*/

	bool operator == (const state& other) const {
		if (pieces.size() != other.pieces.size())
			return false;
		for (int i = 0; i < k; i++) {
			piece* p1 = pieces[i];
			piece* p2 = other.pieces[i];
			if (*p1 != *p2)
				return false;
		}
		return true;
	}

	bool operator != (const state& other) const {
		return !(*this == other);
	}

	bool operator <(const state& other) const {
		return false;
	}

	bool operator > (const state& other) const {
		return false;
	}

	piece* & operator [] (int idx) {
		return pieces[idx];
	}

	piece* const & operator [] (int idx) const {
		return pieces[idx];
	}

	size_t size() {
		return pieces.size();
	}

	vector<vector<state>> generateNextStates() const {
		vector<vector<state>> res;
		res.assign(pieces.size(), vector<state>());
		for (const auto& p : pieces) {
			vector<piece*> ve = p->nextMoves();
			for (const auto& np : ve) {
				res[p->color].push_back(state(*this));
				res[p->color][res[p->color].size() - 1][p->color] = np;
			}
		}
		return res;
	}

	void build(state* prvState = nullptr) {
		if (!prvState)
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++) {
					if (grid[i][j] == WATER || grid[i][j] == BLOCK)
						continue;
					grid[i][j] = WATER;
					color[i][j] = -1;
				}
		else
			for (const auto& p : prvState->pieces)
				for (const auto& init_c : cells[p->color]) {
					cell c = init_c + p->shifts;
					grid[c.first][c.second] = WATER, color[c.first][c.second] = -1;
				}

		for (const auto& p : pieces)
			for (const auto& init_c : cells[p->color]) {
				cell c = init_c + p->shifts;
				grid[c.first][c.second] = p->getChar(), color[c.first][c.second] = p->color;
			}
	}

	//same as dfs/bfs
	bool isFinal() {
		memset(vis, 0, sizeof(bool) * k);
		for (int i = 0; i < m; i++)
			if (color[0][i] >= 0 && !vis[color[0][i]] && pieces[color[0][i]]->getChar() == bridge::CHAR) {
				if (dfs(color[0][i]))
					return true;
			}
		return false;
	}
		
	int heuristic() {
		//return 0;
		memset(vis, 0, sizeof(bool) * k);
		for (int i = 0; i < m; i++)
			if (color[0][i] >= 0 && !vis[color[0][i]] && pieces[color[0][i]]->getChar() == bridge::CHAR)
				Dfs(color[0][i]);

		vector<cell> up_map;
		for (int i = 0; i < k; i++)
			if (vis[i]){
				for (auto c : cells[i])
					up_map.push_back(c + pieces[i]->shifts);
				vis[i] = false;
			}

		for (int i = 0; i < m; i++)
			if (color[n - 1][i] >= 0 && !vis[color[n - 1][i]] && pieces[color[n - 1][i]]->getChar() == bridge::CHAR)
				Dfs(color[n - 1][i]);

		vector<cell> down_map;
		for (int i = 0; i < k; i++)
			if (vis[i])
				for (auto c : cells[i])
					down_map.push_back(c + pieces[i]->shifts);

		int mn = 1e9;
		for (const auto& c1 : up_map)
			for (const auto& c2 : down_map)
				mn = min(mn, cell::manhaten(c1, c2));

		return mn;
	}

	void push_back(piece* p) {
		pieces.push_back(p);
	}

	friend struct hash<state>;
private:
	vector<piece*> pieces;
	//O(K * Log(S) + K * B) with an upper bound of O(N * M)
	bool dfs(int p) {
		vis[p] = true;
		auto it = lower_bound(cells[p].begin(), cells[p].end(), cell(n - 1 - pieces[p]->shifts.first, -1));
		if (it != cells[p].end())
			if ((it->first + pieces[p]->shifts.first) == n - 1)
				return true;

		for (int dir = 0; dir < 4; dir++)
			for (auto idx : bounds[dir][p]) {
				cell nxt = cells[p][idx].dirNeighbour(dir) + pieces[p]->shifts;
				if (cell::valid(nxt, n, m))
					if (grid[nxt.first][nxt.second] == pieces[p]->getChar())
						if (!vis[color[nxt.first][nxt.second]])
							if (dfs(color[nxt.first][nxt.second]))
								return true;
			}

		return false;
	}

	int Dfs(int p) {
		vis[p] = true;
		int sum = 1;
		for (int dir = 0; dir < 4; dir++)
			for (auto idx : bounds[dir][p]) {
				cell nxt = cells[p][idx].dirNeighbour(dir) + pieces[p]->shifts;
				if (cell::valid(nxt, n, m))
					if (grid[nxt.first][nxt.second] == pieces[p]->getChar())
						if (!vis[color[nxt.first][nxt.second]])
							sum += Dfs(color[nxt.first][nxt.second]);
			}
		return sum;
	}

	//O(K * Log(S) + K * B) with an upper bound of O(N * M)
	bool bfs() {
		memset(vis, 0, sizeof(bool) * k);
		queue<int> qu;
		for (int i = 0; i < m; i++)
			if (color[0][i] >= 0 && !vis[color[0][i]] && pieces[color[0][i]]->getChar() == bridge::CHAR)
				qu.push(color[0][i]), vis[color[0][i]] = true;

		while (!qu.empty()) {
			int p = qu.front();
			qu.pop();

			auto it = lower_bound(cells[p].begin(), cells[p].end(), cell(n - 1 - pieces[p]->shifts.first, -1));
			if (it != cells[p].end())
				if ((it->first + pieces[p]->shifts.first) == n - 1)
					return true;

			for (int dir = 0; dir < 4; dir++)
				for (auto idx : bounds[dir][p]) {
					cell nxt = cells[p][idx].dirNeighbour(dir) + pieces[p]->shifts;
					if (cell::valid(nxt, n, m))
						if (grid[nxt.first][nxt.second] == pieces[p]->getChar()) {
							int q = color[nxt.first][nxt.second];
							if (!vis[q]) {
								vis[q] = true;
								qu.push(q);
							}
						}
				}
		}

		return false;
	}
	
};

namespace std {
	template <>
	struct hash<state> {
		size_t operator()(const state& s) const {
			hash<piece*> hasher;
			size_t ret = 0;
			for (const auto& p : s.pieces)
				ret = (ret ^ (hasher(p) << 1)) >> 1;
			return ret;
			size_t res = 1;
			const int prime = 31;
			for (const auto& p : s.pieces)
				res = res * prime + hasher(p);
			return res;
		}
	};
}
#pragma endregion

bool isBounds(const cell& par, const cell& child) {
	return (!cell::valid(child, n, m) && cell::manhaten(par, child) == 1)
		|| color[par.first][par.second] != color[child.first][child.second];
}

//O (N * M)
void split() {
	for (int i = 0; i < k; i++)
		for (int j = 0; j < cells[i].size(); j++) {
			cell& c = cells[i][j];
			vector<cell> ve = c.neighbours();
			for (int dir = 0; dir < 4; dir++)
				if (isBounds(c, ve[dir]))
					bounds[dir][i].push_back(j);
		}
}

state readData(string path) {
	state currentState;
	filebuf fb;
	fb.open(path, ios::in);
	istream is(&fb);

	is >> n >> m;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			is >> grid[i][j];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; k = max(k, color[i][j++]--))
			is >> color[i][j];

	cells.assign(k, vector<cell>());
	for (int i = 0; i < 4; i++)
		bounds[i].assign(cells.size(), vector<int>());

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			if (color[i][j] >= 0)
				cells[color[i][j]].push_back(cell(i, j));
	split();

	for (int i = 0; i < k; i++) {
		cell& c = cells[i][0];
		if (grid[c.first][c.second] == bridge::CHAR)
			currentState.push_back(new bridge(i));
		else
			if (grid[c.first][c.second] == grass::CHAR)
				currentState.push_back(new grass(i));
	}

	return currentState;
}

#pragma endregion

#pragma region gameplay functions
void printGrid() {
	HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	for (int i = 0; i < n; i++, cout << '\n')
		for (int j = 0; j < m; cout << color[i][j++] + 1 << ' ')
			switch (grid[i][j]) {
			case WATER: {
				SetConsoleTextAttribute(hConsole, 3);
				break;
			}
			case BLOCK: {
				SetConsoleTextAttribute(hConsole, 8);
				break;
			}
			case grass::CHAR: {
				SetConsoleTextAttribute(hConsole, 2);
				break;
			}
			case bridge::CHAR: {
				SetConsoleTextAttribute(hConsole, 12);
				break;
			}
			}
	SetConsoleTextAttribute(hConsole, 7);
	for (int j = 0; j < 2 * m; j++)
		cout << '-';
	cout << '\n';
}

//O(S)
void updateGrid(piece* p, const pii& newShifts) {
	for (const auto& init_c : cells[p->color]) {
		cell cc = init_c + p->shifts;
		grid[cc.first][cc.second] = WATER;
		color[cc.first][cc.second] = -1;
	}
	for (const auto& init_c : cells[p->color]) {
		cell nc = init_c + newShifts;
		grid[nc.first][nc.second] = p->getChar();
		color[nc.first][nc.second] = p->color;
	}
}

//O(S)
void updateGridAndPush(state& currentState, piece* p, const pii& newShifts) {
	updateGrid(p, newShifts);
	if (p->getChar() == grass::CHAR)
		currentState[p->color] = new grass(p->color, newShifts);
	if (p->getChar() == bridge::CHAR)
		currentState[p->color] = new bridge(p->color, newShifts);
}

string getMove(const pii& oldShifts, const pii& newShifts) {
	int dx = newShifts.first - oldShifts.first;
	int dy = newShifts.second - oldShifts.second;
	if (dx == 0)
		if (dy > 0)
			return "RIGHT";
		else
			return "LEFT";
	else
		if (dx > 0)
			return "DOWN";
		else
			return "UP";
	return "NULL";
}

int humanPlay() {
	state currentsecondtate = readData("input.txt");
	while (true) {
		::system("cls");
		printGrid();
		cout << "Enter 1 to check your solution or 2 to move a piece:\n";
		int c; cin >> c;
		if (c == 1) {
			if (currentsecondtate.isFinal())
				return cout << "secondolved!!\n", 0;
			cout << "The board is not solved yet!!\n";
		}
		else
			if (c == 2) {
				cout << "Choose a piece, direction and amount:\n";
				int col; cin >> col;
				if (col > k || col < 0) {
					cout << "Wrong input!!\nTry again.\n";
					continue;
				}
				col--;
				debug(DOWN);
				debug(RIGHT);
				debug(UP);
				debug(LEFT);
				int dir; cin >> dir;
				int amount; cin >> amount;
				if (amount < 0) {
					cout << "Wrong input!!\nTry again.\n";
					continue;
				}
				if (dir == DOWN || dir == RIGHT || dir == LEFT || dir == UP)
					if (currentsecondtate[col]->validMove(dir, amount) == amount)
						updateGridAndPush(currentsecondtate, currentsecondtate[col],
							pii(currentsecondtate[col]->shifts.first + ar[dir] * amount,
								currentsecondtate[col]->shifts.second + ac[dir] * amount));
					else
						cout << "Not a valid move!!\nTry again.\n";
				else
					cout << "Wrong input!!\nTry again.\n";
			}
			else
				cout << "Wrong input!!\nTry again.\n";
		::system("pause");
	}

}

vector<pair<int, string>> path;
unordered_map<state, int> dist;

bool dfs(state& s) {
	if (s.isFinal())
		return true;
	dist[s] = true;
	vector<vector<state>> nxt = s.generateNextStates();
	for (int i = 0; i < s.size(); i++) {
		vector<state>& ve = nxt[i];
		for (auto& ns : ve) {
			if (!dist[ns]) {
				updateGrid(s[i], ns[i]->shifts);
				if (dfs(ns))
					return true;
				updateGrid(ns[i], s[i]->shifts);
			}
		}
	}
	return false;
}

bool bfs(state& s) {
	queue<state> qu;
	qu.push(s);
	dist[s] = true;
	while (!qu.empty()) {
		state currentState = qu.front();
		qu.pop();
		currentState.build();
		if (currentState.isFinal())
			return true;
		vector<vector<state>> nxt = currentState.generateNextStates();
		for (auto& ve : nxt) {
			for (auto& ns : ve) {
				if (!dist[ns]){
					dist[ns] = true;
					qu.push(ns);
				}
			}
		}
	}
}

bool Dfs(state& s) {
	stack<state> qu;
	qu.push(s);
	while (!qu.empty()) {
		state currentState = qu.top();
		qu.pop();

		if (dist[currentState])
			continue;
		dist[currentState] = true;

		currentState.build();
		if (currentState.isFinal())
			return true;

		vector<vector<state>> nxt = currentState.generateNextStates();
		for (auto& ve : nxt)
			for (auto& ns : ve)
				if (!dist[ns])
					qu.push(ns);
	}
}

typedef pair<int, state> node;

bool Ucs(state& s) {
	priority_queue<node> pq;
	pq.push(node(0, s));
	dist[s] = 0;
	while (!pq.empty()) {
		node currentNode = pq.top();
		pq.pop();
		state currentState = currentNode.second;
		int currentDist = dist[currentState];
		currentState.build();
		
		if (dist.find(currentState) != dist.end())
			if (currentDist < -currentNode.first)
				continue;
		
		if (currentState.isFinal())
			return true;
		
		vector<vector<state>> nxt = currentState.generateNextStates();
		for (auto& ve : nxt) {
			for (auto& ns : ve) {
				if (dist.find(ns) == dist.end()) {
					dist[ns] = currentDist + 1;
					pq.push(node(-dist[ns], ns));
				}
				else {
					int& Dist = dist[ns];
					if (Dist > currentDist + 1) {
						Dist = currentDist + 1;
						pq.push(node(-Dist, ns));
					}
				}
			}
		}
	}
	return false;
}

unordered_map<state, pii> starDist;

bool Astar(state& s) {
	priority_queue<node> pq;
	pq.push(node(0, s));
	starDist[s] = pii(0, 0);
	while (!pq.empty()) {
		node currentNode = pq.top();
		pq.pop();
		state currentState = currentNode.second;
		int currentG = starDist[currentState].first;
		int currentF = starDist[currentState].second;
		int currentH = currentF - currentG;
		currentState.build();

		if (dist.find(currentState) != dist.end())
			if (currentF < -currentNode.first)
				continue;

		if (currentState.isFinal())
			return true;

		vector<vector<state>> nxt = currentState.generateNextStates();
		for (auto& ve : nxt) {
			for (auto& ns : ve) {
				if (dist.find(ns) == dist.end()) {
					starDist[ns].first = currentG + 1;
					starDist[ns].second = starDist[ns].first + ns.heuristic();
					pq.push(node(-starDist[ns].second, ns));
				}
				else {
					int& lastG = starDist[ns].first;
					int& lastF = starDist[ns].second;
					int newG = currentG + 1;
					int newF = newG + ns.heuristic();
					if (newF > lastF) {
						lastF = newF;
						lastG = newG;
						pq.push(node(-newF, ns));
					}
				}
			}
		}
	}
	return false;
}
#pragma endregion

int main() {
	ios_base::sync_with_stdio(0); cin.tie();

	state currentState = readData("input.txt");
	printGrid();
	
	cout << "solved? " << (Ucs(currentState) ? "Yes!!" : "No!!") << '\n';
	cout << "searched in " << starDist.size() << " states!" << '\n';
	printGrid();
	return 0;
}
