#include <bits/stdc++.h>
using namespace std;

double simp_rate, dist_eps;
int face_num, vertex_num;

int *vertex_id;

bool *bound;

int *pre, *id;
int find(int x) {
	return (pre[x] == x) ? x : pre[x] = find(pre[x]);
}

class vertex {
public:
	double x, y, z;
	double **Q;
	set<int> face_list;
	vertex(double _x = 0., double _y = 0., double _z = 0.): x(_x), y(_y), z(_z) {
		Q = new double*[4];
		for (int i = 0; i < 4; i++) {
			Q[i] = new double[4];
			Q[i][0] = Q[i][1] = Q[i][2] = Q[i][3] = 0;
		}
	}
	void computeQ();
};
vector<vertex> vertices;

class face {
public:
	int v1, v2, v3;
	double a, b, c, d;
	face(int _v1 = 0, int _v2 = 0, int _v3 = 0): v1(_v1), v2(_v2), v3(_v3) {}
	void computePara();
};
vector<face> faces;

vector<set<int> > edges, orig_edges;

inline double dist(int u, int v) {
	double x = vertices[u].x - vertices[v].x, y = vertices[u].y - vertices[v].y, z = vertices[u].z - vertices[v].z;
	return sqrt(x * x + y * y + z * z);
}

void vertex::computeQ() {
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++) Q[i][j] = 0;
	for (auto it = face_list.begin(); it != face_list.end(); it++) {
		double *p = new double[4];
		p[0] = faces[*it].a; p[1] = faces[*it].b; p[2] = faces[*it].c; p[3] = faces[*it].d;
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++) Q[i][j] += p[i] * p[j];
		delete [] p;
	}
}

void face::computePara() {
	assert(v1 == pre[v1]); assert(v2 == pre[v2]); assert(v3 == pre[v3]);
	double x1 = vertices[v2].x - vertices[v1].x, y1 = vertices[v2].y - vertices[v1].y, z1 = vertices[v2].z - vertices[v1].z;
	double x2 = vertices[v3].x - vertices[v1].x, y2 = vertices[v3].y - vertices[v1].y, z2 = vertices[v3].z - vertices[v1].z;
	a = y1 * z2 - y2 * z1; b = z1 * x2 - z2 * x1; c = x1 * y2 - x2 * y1;
	double len = sqrt(a * a + b * b + c * c);
	a /= len; b /= len; c /= len;
	d = -(a * vertices[v1].x + b * vertices[v1].y + c * vertices[v1].z);
	/*if (fabs(a * vertices[v2].x + b * vertices[v2].y + c * vertices[v2].z + d) > 1e-9) {
		fprintf(stderr, "%d %d %d\n", v1, v2, v3);
		fprintf(stderr, "%.10lf %.10lf %.10lf\n", vertices[v1].x, vertices[v1].y, vertices[v1].z);
		fprintf(stderr, "%.10lf %.10lf %.10lf\n", vertices[v2].x, vertices[v2].y, vertices[v2].z);
		fprintf(stderr, "%.10lf %.10lf %.10lf\n", vertices[v3].x, vertices[v3].y, vertices[v3].z);
		fprintf(stderr, "%.10lf\n", len);
	} else if (fabs(a * vertices[v3].x + b * vertices[v3].y + c * vertices[v3].z + d) > 1e-9) {
		fprintf(stderr, "%d %d %d\n", v1, v2, v3);
		fprintf(stderr, "%.10lf %.10lf %.10lf\n", vertices[v1].x, vertices[v1].y, vertices[v1].z);
		fprintf(stderr, "%.10lf %.10lf %.10lf\n", vertices[v2].x, vertices[v2].y, vertices[v2].z);
		fprintf(stderr, "%.10lf %.10lf %.10lf\n", vertices[v3].x, vertices[v3].y, vertices[v3].z);
		fprintf(stderr, "%.10lf\n", len);
	}*/
	//assert(fabs(a * vertices[v2].x + b * vertices[v2].y + c * vertices[v2].z + d) < 1e-9);
	//assert(fabs(a * vertices[v3].x + b * vertices[v3].y + c * vertices[v3].z + d) < 1e-9);
}

class ver_pair {
public:
	int u, v;
	double err;
	double x, y, z;
	int time_cnt;
	ver_pair(int _u, int _v, double _err, double _x, double _y, double _z, int _time): u(_u), v(_v), err(_err), x(_x), y(_y), z(_z), time_cnt(_time) {}
};
inline bool operator <(const ver_pair &a, const ver_pair &b) {
	return a.err > b.err;
}

priority_queue<ver_pair> heap;
map<pair<int, int>, int> timer;

double computeMinCost(int u, int v, double *_x, double *_y, double *_z) {
	assert(pre[u] == u); assert(pre[v] == v);
	double **A = new double*[4];
	for (int i = 0; i < 4; i++) {
		A[i] = new double[4];
		for (int j = 0; j < 4; j++) A[i][j] = vertices[u].Q[i][j] + vertices[v].Q[i][j];
	}
	A[3][0] = A[3][1] = A[3][2] = 0; A[3][3] = 1;
	double *b = new double[4];
	b[0] = b[1] = b[2] = 0; b[3] = 1;
	double *x = new double[4];
	x[0] = x[1] = x[2] = x[3] = 0;
	bool flag = true;
	for (int i = 0; i < 4; i++) {
		int maxLine = i;
		for (int j = i + 1; j < 4; j++)
			if (fabs(A[j][i]) > fabs(A[maxLine][i])) maxLine = j;
		for (int j = i; j < 4; j++) swap(A[i][j], A[maxLine][j]);
		swap(b[i], b[maxLine]);
		double t = A[i][i];
		if (fabs(t) < 1e-10) {
			x[0] = (vertices[u].x + vertices[v].x) / 2.;
			x[1] = (vertices[u].y + vertices[v].y) / 2.;
			x[2] = (vertices[u].z + vertices[v].z) / 2.;
			x[3] = 1.;
			flag = false;
			break;
		}
		for (int j = i; j < 4; j++) A[i][j] /= t;
		b[i] /= t;
		for (int j = i + 1; j < 4; j++) if (fabs(A[j][i]) > 1e-8) {
			t = A[j][i];
			for (int k = i; k < 4; k++) A[j][k] -= A[i][k] * t;
			b[j] -= b[i] * t;
		}
	}
	if (flag) {
		for (int i = 3; i >= 0; i--) {
			x[i] = b[i];
			for (int k = i + 1; k < 4; k++) x[i] -= A[i][k] * x[k];
		}
	}
	assert(fabs(x[3] - 1.) < 1e-8);
	*_x = x[0]; *_y = x[1]; *_z = x[2];
	double cost = 0;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			cost += x[i] * x[j] * (vertices[u].Q[i][j] + vertices[v].Q[i][j]);
		delete [] A[i];
	}
	for (auto it = vertices[u].face_list.begin(); it != vertices[u].face_list.end(); it++) {
		if (faces[*it].v1 == v || faces[*it].v2 == v || faces[*it].v3 == v) continue;
		int v1 = faces[*it].v1, v2 = faces[*it].v2, v3 = faces[*it].v3;
		double x1, y1, z1, x2, y2, z2;
		double na, nb, nc, nna, nnb, nnc;
		if (faces[*it].v1 == u) {
			x1 = vertices[v2].x - vertices[v1].x; y1 = vertices[v2].y - vertices[v1].y; z1 = vertices[v2].z - vertices[v1].z;
			x2 = vertices[v3].x - vertices[v1].x; y2 = vertices[v3].y - vertices[v1].y; z2 = vertices[v3].z - vertices[v1].z;
			na = y1 * z2 - y2 * z1; nb = z1 * x2 - z2 * x1; nc = x1 * y2 - x2 * y1;
			x1 = vertices[v2].x - x[0]; y1 = vertices[v2].y - x[1]; z1 = vertices[v2].z - x[2];
			x2 = vertices[v3].x - x[0]; y2 = vertices[v3].y - x[1]; z2 = vertices[v3].z - x[2];
			nna = y1 * z2 - y2 * z1; nnb = z1 * x2 - z2 * x1; nnc = x1 * y2 - x2 * y1;
		} else if (faces[*it].v2 == u) {
			x1 = vertices[v1].x - vertices[v2].x; y1 = vertices[v1].y - vertices[v2].y; z1 = vertices[v1].z - vertices[v2].z;
			x2 = vertices[v3].x - vertices[v2].x; y2 = vertices[v3].y - vertices[v2].y; z2 = vertices[v3].z - vertices[v2].z;
			na = y1 * z2 - y2 * z1; nb = z1 * x2 - z2 * x1; nc = x1 * y2 - x2 * y1;
			x1 = vertices[v1].x - x[0]; y1 = vertices[v1].y - x[1]; z1 = vertices[v1].z - x[2];
			x2 = vertices[v3].x - x[0]; y2 = vertices[v3].y - x[1]; z2 = vertices[v3].z - x[2];
			nna = y1 * z2 - y2 * z1; nnb = z1 * x2 - z2 * x1; nnc = x1 * y2 - x2 * y1;
		} else {
			assert(faces[*it].v3 == u);
			x1 = vertices[v1].x - vertices[v3].x; y1 = vertices[v1].y - vertices[v3].y; z1 = vertices[v1].z - vertices[v3].z;
			x2 = vertices[v2].x - vertices[v3].x; y2 = vertices[v2].y - vertices[v3].y; z2 = vertices[v2].z - vertices[v3].z;
			na = y1 * z2 - y2 * z1; nb = z1 * x2 - z2 * x1; nc = x1 * y2 - x2 * y1;
			x1 = vertices[v1].x - x[0]; y1 = vertices[v1].y - x[1]; z1 = vertices[v1].z - x[2];
			x2 = vertices[v2].x - x[0]; y2 = vertices[v2].y - x[1]; z2 = vertices[v2].z - x[2];
			nna = y1 * z2 - y2 * z1; nnb = z1 * x2 - z2 * x1; nnc = x1 * y2 - x2 * y1;
		}
		if (na * nna + nb * nnb + nc * nnc <= 0) {
			cost = 1e5;
			//fprintf(stderr, "inverse\n");
			break;
		}
	}
	for (auto it = vertices[v].face_list.begin(); it != vertices[v].face_list.end(); it++) {
		if (faces[*it].v1 == u || faces[*it].v2 == u || faces[*it].v3 == u) continue;
		int v1 = faces[*it].v1, v2 = faces[*it].v2, v3 = faces[*it].v3;
		double x1, y1, z1, x2, y2, z2;
		double na, nb, nc, nna, nnb, nnc;
		if (faces[*it].v1 == v) {
			x1 = vertices[v2].x - vertices[v1].x; y1 = vertices[v2].y - vertices[v1].y; z1 = vertices[v2].z - vertices[v1].z;
			x2 = vertices[v3].x - vertices[v1].x; y2 = vertices[v3].y - vertices[v1].y; z2 = vertices[v3].z - vertices[v1].z;
			na = y1 * z2 - y2 * z1; nb = z1 * x2 - z2 * x1; nc = x1 * y2 - x2 * y1;
			x1 = vertices[v2].x - x[0]; y1 = vertices[v2].y - x[1]; z1 = vertices[v2].z - x[2];
			x2 = vertices[v3].x - x[0]; y2 = vertices[v3].y - x[1]; z2 = vertices[v3].z - x[2];
			nna = y1 * z2 - y2 * z1; nnb = z1 * x2 - z2 * x1; nnc = x1 * y2 - x2 * y1;
		} else if (faces[*it].v2 == v) {
			x1 = vertices[v1].x - vertices[v2].x; y1 = vertices[v1].y - vertices[v2].y; z1 = vertices[v1].z - vertices[v2].z;
			x2 = vertices[v3].x - vertices[v2].x; y2 = vertices[v3].y - vertices[v2].y; z2 = vertices[v3].z - vertices[v2].z;
			na = y1 * z2 - y2 * z1; nb = z1 * x2 - z2 * x1; nc = x1 * y2 - x2 * y1;
			x1 = vertices[v1].x - x[0]; y1 = vertices[v1].y - x[1]; z1 = vertices[v1].z - x[2];
			x2 = vertices[v3].x - x[0]; y2 = vertices[v3].y - x[1]; z2 = vertices[v3].z - x[2];
			nna = y1 * z2 - y2 * z1; nnb = z1 * x2 - z2 * x1; nnc = x1 * y2 - x2 * y1;
		} else {
			assert(faces[*it].v3 == v);
			x1 = vertices[v1].x - vertices[v3].x; y1 = vertices[v1].y - vertices[v3].y; z1 = vertices[v1].z - vertices[v3].z;
			x2 = vertices[v2].x - vertices[v3].x; y2 = vertices[v2].y - vertices[v3].y; z2 = vertices[v2].z - vertices[v3].z;
			na = y1 * z2 - y2 * z1; nb = z1 * x2 - z2 * x1; nc = x1 * y2 - x2 * y1;
			x1 = vertices[v1].x - x[0]; y1 = vertices[v1].y - x[1]; z1 = vertices[v1].z - x[2];
			x2 = vertices[v2].x - x[0]; y2 = vertices[v2].y - x[1]; z2 = vertices[v2].z - x[2];
			nna = y1 * z2 - y2 * z1; nnb = z1 * x2 - z2 * x1; nnc = x1 * y2 - x2 * y1;
		}
		if (na * nna + nb * nnb + nc * nnc <= 0) {
			cost = 1e5;
			//fprintf(stderr, "inverse\n");
			break;
		}
	}
	delete [] A;
	delete [] x;
	delete [] b;
	return cost;
}

inline bool valid(const ver_pair &edge) {
	if (edge.u != pre[edge.u]) return false;
	if (edge.v != pre[edge.v]) return false;
	if (bound[edge.u] || bound[edge.v]) return false;
	if (edge.err > 1e4) return false;
	return edge.time_cnt == timer[make_pair(edge.u, edge.v)];
}

void simplify() {
	int u, v;
	for (u = 0; u < vertex_num; u++)
		for (auto it = edges[u].begin(); it != edges[u].end(); it++) {
			v = *it;
			if (u > v) continue;
			timer[make_pair(u, v)] = 1;
			double x, y, z;
			double err = computeMinCost(u, v, &x, &y, &z);
			heap.push(ver_pair(u, v, err, x, y, z, 1));
		}
	
	int target = (int)((1 - simp_rate) * face_num) >> 1;
	//fprintf(stderr, "%d\n", target);
	for (int cnt = 0; cnt < target; cnt++) {
		//fprintf(stderr, "%d\n", cnt);
		while (!heap.empty() && !valid(heap.top())) heap.pop();
		if (heap.empty()) {
			fprintf(stderr, "%d\n", cnt);
			break;
		}
		ver_pair edge = heap.top(); heap.pop();
		u = edge.u; v = edge.v;
		//if (cnt < 10) fprintf(stderr, "%d %d\n", u, v);
		assert(u < v);
		pre[v] = u;
		vertices[u].x = edge.x; vertices[u].y = edge.y; vertices[u].z = edge.z;
		edges[u].erase(v); orig_edges[u].erase(v);
		for (auto it = vertices[u].face_list.begin(); it != vertices[u].face_list.end();) {
			if (faces[*it].v1 == v || faces[*it].v2 == v || faces[*it].v3 == v) {
				int f = *it;
				it++;
				vertices[faces[f].v1].face_list.erase(f);
				vertices[faces[f].v2].face_list.erase(f);
				vertices[faces[f].v3].face_list.erase(f);
			} else {
				faces[*it].computePara();
				it++;
			}
		}
		for (auto it = vertices[v].face_list.begin(); it != vertices[v].face_list.end(); it++) {
			if (faces[*it].v1 != u && faces[*it].v2 != u && faces[*it].v3 != u) {
				if (faces[*it].v1 == v) {
					faces[*it].v1 = u;
					edges[u].insert(faces[*it].v2); orig_edges[u].insert(faces[*it].v2);
					edges[faces[*it].v2].erase(v); orig_edges[faces[*it].v2].erase(v);
					edges[faces[*it].v2].insert(u); orig_edges[faces[*it].v2].insert(u);
					edges[u].insert(faces[*it].v3); orig_edges[u].insert(faces[*it].v3);
					edges[faces[*it].v3].erase(v); orig_edges[faces[*it].v3].erase(v);
					edges[faces[*it].v3].insert(u); orig_edges[faces[*it].v3].insert(u);
				} else if (faces[*it].v2 == v) {
					faces[*it].v2 = u;
					edges[u].insert(faces[*it].v1); orig_edges[u].insert(faces[*it].v1);
					edges[faces[*it].v1].erase(v); orig_edges[faces[*it].v1].erase(v);
					edges[faces[*it].v1].insert(u); orig_edges[faces[*it].v1].insert(u);
					edges[u].insert(faces[*it].v3); orig_edges[u].insert(faces[*it].v3);
					edges[faces[*it].v3].erase(v); orig_edges[faces[*it].v3].erase(v);
					edges[faces[*it].v3].insert(u); orig_edges[faces[*it].v3].insert(u);
				} else if (faces[*it].v3 == v) {
					faces[*it].v3 = u;
					edges[u].insert(faces[*it].v1); orig_edges[u].insert(faces[*it].v1);
					edges[faces[*it].v1].erase(v); orig_edges[faces[*it].v1].erase(v);
					edges[faces[*it].v1].insert(u); orig_edges[faces[*it].v1].insert(u);
					edges[u].insert(faces[*it].v2); orig_edges[u].insert(faces[*it].v2);
					edges[faces[*it].v2].erase(v); orig_edges[faces[*it].v2].erase(v);
					edges[faces[*it].v2].insert(u); orig_edges[faces[*it].v2].insert(u);
				}
				faces[*it].computePara();
				vertices[u].face_list.insert(*it);
			}
		}
		vertices[u].computeQ();
		for (auto it = edges[u].begin(); it != edges[u].end(); it++) {
			v = *it;
			if (pre[v] != v) {
				edges[u].erase(v);
				continue;
			}
			vertices[v].computeQ();
			assert(u != v);
			int new_time = ++timer[(u < v) ? make_pair(u, v) : make_pair(v, u)];
			if ((orig_edges[u].find(v) != orig_edges[u].end()) || (dist(u, v) < dist_eps)) {
				double x, y, z;
				double err = computeMinCost(u, v, &x, &y, &z);
				heap.push((u < v) ? ver_pair(u, v, err, x, y, z, new_time) : ver_pair(v, u, err, x, y, z, new_time));
			}
		}
	}
}

int D;

inline bool cmp(const int &u, const int &v) {
	return (D == 0) ? vertices[u].x < vertices[v].x : ((D == 1) ? vertices[u].y < vertices[v].y : vertices[u].z < vertices[v].z);
}

void find_close_pairs(int l, int r, int dim) {
	if (dim == 3) dim = 0;
	int mid = (l + r) >> 1;
	D = dim;
	sort(vertex_id + l, vertex_id + r, cmp);
	double pos;
	if (dim == 0)
		pos = vertices[vertex_id[mid]].x;
	else if (dim == 1)
		pos = vertices[vertex_id[mid]].y;
	else
		pos = vertices[vertex_id[mid]].z;
	for (int i = mid; i >= l; i--) {
		if (dim == 0) {
			if (vertices[vertex_id[i]].x + dist_eps < pos) break;
		} else if (dim == 1) {
			if (vertices[vertex_id[i]].y + dist_eps < pos) break;
		} else {
			if (vertices[vertex_id[i]].z + dist_eps < pos) break;
		}
		for (int j = mid + 1; j < r; j++) {
			if (dim == 0) {
				if (vertices[vertex_id[j]].x - dist_eps > pos) break;
			} else if (dim == 1) {
				if (vertices[vertex_id[j]].y - dist_eps > pos) break;
			} else {
				if (vertices[vertex_id[j]].z - dist_eps > pos) break;
			}
			if (dist(vertex_id[i], vertex_id[j]) < dist_eps) {
				edges[vertex_id[i]].insert(vertex_id[j]);
				edges[vertex_id[j]].insert(vertex_id[i]);
			}
		}
	}
	if (l + 1 < mid) find_close_pairs(l, mid, dim + 1);
	if (mid + 1 < r) find_close_pairs(mid, r, dim + 1);
}

int main(int argc, char *argv[]) {
	if (argc < 5) {
		puts("Usage: mesh_simp [input model] [output model] [simp rate] [dist threshold]");
		return 0;
	}
	freopen(argv[1], "r", stdin);
	freopen(argv[2], "w", stdout);
	simp_rate = atof(argv[3]);
	dist_eps = atof(argv[4]);
	ios::sync_with_stdio(false);
	
	char c;
	while ((c = getchar()) && (c != EOF)) {
		if (c == '#') {
			while ((c = getchar()) && (c != EOF) && (c != '\n'));
			continue;
		}
		if (c == 'v') {
			double x, y, z;
			scanf("%lf%lf%lf", &x, &y, &z);
			vertices.push_back(vertex(x, y, z));
			vertex_num++;
		} else if (c == 'f') {
			int u, v, w;
			scanf("%d%d%d", &u, &v, &w);
			u--; v--; w--;
			faces.push_back(face(u, v, w));
			assert((u < vertex_num) && (v < vertex_num) && (w < vertex_num));
			vertices[u].face_list.insert(face_num);
			vertices[v].face_list.insert(face_num);
			vertices[w].face_list.insert(face_num);
			face_num++;
		}
	}
	
	edges.resize(vertex_num);
	orig_edges.resize(vertex_num);
	for (int i = 0; i < face_num; i++) {
		int u = faces[i].v1, v = faces[i].v2, w = faces[i].v3;
		assert((u != v) && (v != w) && (w != u));
		edges[u].insert(v); edges[u].insert(w);
		edges[v].insert(u); edges[v].insert(w);
		edges[w].insert(u); edges[u].insert(v);
		orig_edges[u].insert(v); orig_edges[u].insert(w);
		orig_edges[v].insert(u); orig_edges[v].insert(w);
		orig_edges[w].insert(u); orig_edges[u].insert(v);
	}
	
	pre = new int[vertex_num];
	for (int i = 0; i < vertex_num; i++) pre[i] = i;
	for (int i = 0; i < face_num; i++) faces[i].computePara();
	for (int i = 0; i < vertex_num; i++) vertices[i].computeQ();
	
	bound = new bool[vertex_num];
	for (int i = 0; i < vertex_num; i++) bound[i] = false;
	for (int u = 0; u < vertex_num; u++) if (!bound[u])
		for (auto it = edges[u].begin(); it != edges[u].end(); it++) {
			int v = *it;
			int face_cnt = 0;
			for (auto it = vertices[u].face_list.begin(); it != vertices[u].face_list.end(); it++)
				if (faces[*it].v1 == v || faces[*it].v2 == v || faces[*it].v3 == v) face_cnt++;
			if (face_cnt < 2)
				bound[u] = bound[v] = true;
		}
	
	vertex_id = new int[vertex_num];
	for (int i = 0; i < vertex_num; i++) vertex_id[i] = i;
	find_close_pairs(0, vertex_num, 0);
	
	simplify();
	
	id = new int[vertex_num];
	int cnt = 0;
	for (int i = 0; i < vertex_num; i++)
		if (pre[i] == i) {
			id[i] = cnt++;
			printf("v %.4lf %.4lf %.4lf\n", vertices[i].x, vertices[i].y, vertices[i].z);
		} else id[i] = -1;
	for (int i = 0; i < face_num; i++)
		if (pre[faces[i].v1] == faces[i].v1 && pre[faces[i].v2] == faces[i].v2 && pre[faces[i].v3] == faces[i].v3)
			printf("f %d %d %d\n", id[faces[i].v1] + 1, id[faces[i].v2] + 1, id[faces[i].v3] + 1);
	
	delete [] vertex_id;
	delete [] pre;
	delete [] id;
	delete [] bound;
	
	return 0;
}
