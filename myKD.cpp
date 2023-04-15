/*

PROJECT: KD-TREE
AUTHOR: ALEXANDER SEBASTIÁN GÓMEZ DEL CARPIO
DATE: 04/28/2022

[EN]
In this project we can see a simple implementation of the KD-Tree with templates for the data type and for the number of dimensions, it makes use
of a "Node" structure and another one of the KDTree with their respective functions. We also find global functions of a structure like sortV() that sorts
a vector and armarTree() that shapes the tree with initial values.

[ES]
En el presente código podemos ver un implementación simple del KD-Tree con templates para el tipo de dato y para la cantidad de dimensiones, este hace uso
de una estructura "Node" y otra del mismo KDTree con sus respectivas funciones. Además encontramos funciones libres de una estructura como sortV() que ordena 
un vector y armarTree() que con valores de inicio se da forma al árbol.
*/

#include<iostream>
#include <algorithm> 
#include<vector>
#include<cmath>
#include<math.h>

using namespace std;

template<class T, int dim>
struct Node {			//Node structure
	Node<T, dim>* right = NULL;			//Pointers to the sides
	Node<T, dim>* left = NULL;
	vector<T> val;
	int nivel = 0;			//Its level in the tree
	Node() {}
	~Node() {}
};

template<class T, int dim>		//Template T: Datatype and dim: Number of dimensions
struct KDTree {				
	Node<T, dim>* root = NULL;
	KDTree() {}
	~KDTree() {}

	void DFS(Node<T, dim>*& ptr) {		//Prints the tree in any dimension
		if (ptr == nullptr) {
			return;
		}
		cout << endl << ptr->val[0] << " " << ptr->val[1] << " Nivel: " << ptr->nivel;
		DFS(ptr->right);
		DFS(ptr->left);
	}

	double distance(Node<T,dim>* a,vector<T> b) {		//Discover the distance between the node and coordinates
		double dist = 0;
		double suma = 0;
		for (int i = 0; i < dim; i++)
		{
			suma += pow(a->val[i] - b[i], 2);
		}
		dist = sqrt(suma);
		return dist;
	}

	void rangeQuery2(vector<T> query_i,vector<T> query_f,Node<T,dim>* ptr, vector<Node<T,dim>*> &respuesta) {		//Recursive function to find put all the tree points inside a range
		if (ptr) {
			int discriminador = ptr->nivel % dim;
			if (query_i[discriminador] <= ptr->val[discriminador] && ptr->val[discriminador] <= query_f[discriminador]) {		//If the coordinate is inside the range, we keep searching to the left as well as to the right
				bool contiene = true;
				for (int i = 0; contiene && i < dim; i++) {							//Check if all the dimensions are inside the range
					if (query_i[i] > ptr->val[i] || ptr->val[i] > query_f[i]) {
						contiene = false;
					}
				}
				if (contiene) {
					respuesta.push_back(ptr);
				}
				rangeQuery2(query_i, query_f, ptr->left, respuesta);
				rangeQuery2(query_i, query_f, ptr->right, respuesta);
			}
			else if (query_i[discriminador] > ptr->val[discriminador]) {		//If it's not in the range, we look for a point that is
				rangeQuery2(query_i, query_f, ptr->right, respuesta);
			}
			else {
				rangeQuery2(query_i, query_f, ptr->left, respuesta);
			}
		}
		return;
	}

	void rangeQuery() {		//It gets the search range and then it calls 'rangeQuery2()' in order to return the points in that interval.
		vector<T> query_i;
		vector<T> query_f;
		vector<Node<T,dim>*> respuesta;
		T v1,v2;
		cout << "\nIngresar intervalos (de menor a mayor) de busqueda\n";  //COUT the interval
		for (int i = 0; i < dim; i++)
		{
			cout << "\tPara dimension "<< i << ": ";
			cin >> v1 >> v2;
			query_i.push_back(v1); //Lower limit
			query_f.push_back(v2); //Upper limit
		}
		rangeQuery2(query_i,query_f,root,respuesta);
		cout << "Los nodos en ese rango son: ";
		for (int i = 0; i < respuesta.size(); i++){
			cout << "\n\t";
			for (int j = 0; j < dim; j++) {
				cout<<respuesta[i]->val[j]<<" ";
			}
		}
	}

	void NN2(vector<T> pto,Node<T,dim>* now,Node<T,dim>*& candidate, int depth, double& bestD=-1) {			//Find the closest point to the point entered in NN with recursion.
		if (!now) {
			return;
		}
		if (bestD==-1 || distance(now, pto) < bestD) {			//If the current distance is better than the one already saved, we replace it
			bestD = distance(now, pto);
			candidate = now;
		}

		int axis = depth % dim;
		bool right = false;
		if (pto[axis] < now->val[axis]) {				//We see if the point is on the right or left side
			right = false;
			NN2(pto,now->left,candidate,depth++,bestD);
		}
		else {
			right = true;
			NN2(pto,now->right, candidate, depth++, bestD);
		}
		if (fabs(now->val[axis] - pto[axis]) < candidate->val[axis]) {			//If the current distance making a circle exceeds the side to be ignored
			if (right) {
				NN2(pto,now->left, candidate, depth++, bestD);
			}
			else {
				NN2(pto,now->right, candidate, depth++,bestD);
			}
		}
	}

	void NN() {		//Gets a non-existent point in the tree and calls NN().
		vector<T> query;
		T v;
		cout << "\nIngrese las coordenadas para buscar su NN: ";
		for (int i = 0; i < dim; i++)
		{
			cin >> v;
			query.push_back(v);
		}
		Node<T, dim>* winner=nullptr;
		double best_distance = -1;
		NN2(query,root,winner,0,best_distance);
		cout << "Su nodo mas cercano es: ";
		for (int i = 0; winner && i < dim; i++)
		{
			cout << winner->val[i]<<" ";
		}

		return;
	}
	void KNN2(vector<T> query, vector<pair<double,Node<T, dim>*>> &ans, Node<T, dim>* ptr) {		// It recursively traverses the tree and finds the K closest points to the point entered in KNN.
		if (!ptr) {
			return;
		}
		ans.push_back(make_pair(distance(ptr,query),ptr));
		KNN2(query,ans,ptr->right);		
		KNN2(query,ans,ptr->left);
	}

	void KNN() {		// Gets a non-existent point in the tree and calls KNN().
		vector<T> query;
		vector<pair<double,Node<T,dim>*>> ans;
		T v;
		int K;
		cout << "\n\nIngrese numero de vecinos a encontrar: ";
		cin >> K;
		cout << "\nIngrese las coordenadas para buscar sus KNN: ";
		for (int i = 0; i < dim; i++)
		{
			cin >> v;
			query.push_back(v);
		}

		KNN2(query, ans, root);

		sort(ans.begin(),ans.end());		//Sorts all distances to the entered point
		
		for (int i = 0; i < ans.size() && i < K; i++) {
			cout << "\n\t";
			for (int j = 0; j < dim; j++) {
				cout<<ans[i].second->val[j]<<" ";
			}
		}

		return;
	}

	Node<T, dim>* search(Node<T, dim>* root2, Node<T, dim>*& preptr, vector<T>& query, int n) {			//Search a node in the tree
		Node<T, dim>* resultado = nullptr;

		if (query.size() == 0) {
			cout << "\nInserta las coordenadas actuales: ";
			T v;
			for (int i = 0; i < dim; i++) {
				cin >> v;
				query.push_back(v);
			}
		}

		int nivel = n; 
		Node<T, dim>* ptr = root2; 
		Node<T, dim>* tmp;   

		while (ptr) {
			if (ptr->val[nivel % dim] > query[nivel % dim]) {   // To the left if the value to look for is less
				preptr = ptr;
				ptr = ptr->left;
			}
			else if (ptr->val[nivel % dim] < query[nivel % dim]) {   //To the right if it is greater
				preptr = ptr;
				ptr = ptr->right;
			}
			else {						//if they are the same
				resultado = ptr;
				for (int i = 0; resultado && i < dim; i++) {		//We see if it is the point that is being sought
					if (ptr->val[i] != query[i]) {
						resultado = nullptr;
					}
				}
				if (resultado) return resultado;
				else {												//If it is not, we go both ways.
					tmp = ptr;
					resultado = search(ptr->left, tmp, query, nivel + 1);
					if (resultado) {
						preptr = tmp;
						return resultado;
					}
					resultado = search(ptr->right, tmp, query, nivel + 1);
					if (resultado) {
						preptr = tmp;
						return resultado;
					}
				}
			}
			nivel++;
		}
		return resultado;
	}

	bool insert() {				//Inserts a node in the tree
		vector<T> query;
		cout << "\nInserta las nuevas coordenadas: ";
		T v;
		for (int i = 0; i < dim; i++) {
			cin >> v;
			query.push_back(v);
		}
		Node<T, dim>* preptr = root;						//It will point to the node one level above the insertion location of the new node
		Node<T, dim>* ptr = search(root, preptr, query, 0);			//If the node is already in the tree we no longer insert it

		Node<T, dim>* nuevo = new Node<T, dim>;
		if (ptr) {
			cout << "El valor ya existe, ingrese otro nodo\n";
			return false;
		}
		else {
			nuevo->val = query;
			if (!preptr) {				//If the tree is empty
				root = nuevo;
			}
			else {
				nuevo->nivel = preptr->nivel + 1;
				if (query[preptr->nivel % dim] > preptr->val[preptr->nivel % dim]) {		//We see if we insert it to the right or to the left
					preptr->right = nuevo;
				}
				else {
					preptr->left = nuevo;
				}
			}
			return true;
		}
	}



	void findMin(Node<T, dim>*& min, Node<T, dim>* tmp, int d, Node<T, dim>*& preptr, Node<T, dim>* pretmp) {			//Find the minimum value in a specific dimension in a subtree
		if (!tmp) {
			return;
		}
		if (tmp->val[d] <= min->val[d]) {
			preptr = pretmp;
			min = tmp;
		}
		findMin(min, tmp->right, d, preptr, tmp);
		findMin(min, tmp->left, d, preptr, tmp);

	}

	void borrar2(Node<T, dim>* ptr, Node<T, dim>* preptr, Node<T, dim>* min = nullptr) {	//Recursive function to delete a node
		if (ptr) {	//It is found
			if (ptr->right) {				//If it is not a leaf and it has a child on the right
				min = ptr->right;
				preptr = ptr;
				findMin(min, ptr->right, ptr->nivel % dim, preptr, ptr);		//We must find the minimum in the discriminant dimension and in its subtree
				ptr->val = min->val;
				borrar2(min, preptr);
			}
			else if (ptr->left) {			//He has a child on the left who will later be the children on the right
				min = ptr->left;
				preptr = ptr;
				findMin(min, ptr->left, ptr->nivel % dim, preptr, ptr);
				ptr->val = min->val;
				ptr->right = ptr->left;
				ptr->left = nullptr;
				borrar2(min, preptr);
			}
			else {							//If it is a leaf, it is only deleted
				if (preptr->right == ptr) {
					preptr->right = nullptr;
				}
				else {
					preptr->left = nullptr;
				}
			}
		}
		return;
	}
	bool borrar() {				//Gets the value to delete
		vector<T> query;
		cout << "\nInserta las coordenadas del nodo a borrar: ";
		T v;
		for (int i = 0; i < dim; i++) {
			cin >> v;
			query.push_back(v);
		}
		Node<T, dim>* preptr = root;
		Node<T, dim>* ptr = search(root, preptr, query, 0);
		Node<T, dim>* min;

		if (ptr) {
			borrar2(ptr, preptr);
			return true;
		}
		cout << "No existe el nodo\n";
		return false;
	}

	void update() {			//Updates a value in the tree
		vector<T> nuevas;
		vector<T> query;
		if (borrar()) {
			while(!insert()){}
		}
		return;
	}
};


template<class T, int dim>
struct compare						//Helps the sort function of a vector to compare and order it
{
	int d = 0; //Dimension to compare
	compare(int _d) :d(_d) {}
	inline bool operator() (const Node<T, dim>* struct1, const Node<T, dim>* struct2)
	{
		return (struct1->val[d] < struct2->val[d]);
	}
};


template<class T, int dim>
void sortV(vector<Node<T, dim>*>& v, int d) {   //Sorts in a specific dimension
	std::sort(v.begin(), v.end(), compare<T, dim>(d));
	return;
}

template<class T, int dim>
void armarTree(KDTree<T, dim>& KD, vector<Node<T, dim>*>& v, Node<T, dim>*& ptr, int d) {		//With data entered at the start, the start tree is built
	int size = v.size();
	int i_nodo = size / 2;
	if (size == 0) {
		return;
	}
	sortV(v, d % dim);		//Order the vector of nodes in the corresponding dimension
	ptr = v[i_nodo];
	v[i_nodo]->nivel = d;
	vector<Node<T, dim>*> v_left;
	copy(v.begin(), v.begin() + i_nodo, back_inserter(v_left));		//We cut the vector in two to call armarTree() for its respective construction
	vector<Node<T, dim>*> v_right;
	copy(v.begin() + i_nodo + 1, v.end(), back_inserter(v_right));

	armarTree(KD, v_left, ptr->left, d + 1);
	armarTree(KD, v_right, ptr->right, d + 1);

	return;
}


int main() {
	const int n_nodos = 8;	//No of nodes
	const int d = 2;	//dimension
	vector<Node<int, d>*> v_inicio;			//Stores the startup nodes
	Node<int, d>* nodos = new Node<int, d>[n_nodos]; //Creating startup nodes
	KDTree<int, d> KDT;

	//INSTANCE
	int instancia[n_nodos][d] = {{7,7},{1,2},{3,4},{6,9},{9,7},{18,10},{7,2},{1,10}};
	for (int i = 0; i < n_nodos; i++) {
		v_inicio.push_back(&nodos[i]);
		for (int j = 0; j < d; j++) {
			v_inicio[i]->val.push_back(instancia[i][j]);
		}
	}

	//LOOPS TO REQUEST VALUES FROM THE USER IN THE TERMINAL//////////////////////////////////////////////
	/*int v;
	for (int i = 0; i < n_nodos; i++) {
		v_inicio.push_back(&nodos[i]);	   //Filling vector with initial coordinates
		for (int j = 0; j < d; j++) {		
			cin >> v;
			v_inicio[i]->val.push_back(v);
		}
	}*/
	//////////////////////////////////////////////////////////////////////////////////////

	armarTree(KDT, v_inicio, KDT.root, 0);		//we put together the tree

	KDT.DFS(KDT.root);			//We print the initial tree

	int tests = 1;			//Number of tests
	for(int i=0;i<tests ;i++){	//UNCOMMENT TO TEST EACH FUNCTION----------------

		KDT.KNN();				//The K nearest neighbors	

		//KDT.NN();						//the nearest neighbor
		
		
		//KDT.rangeQuery();				 //Points within a range
		

		/*KDT.insert();					//Delete and insert functions have the search function inside
		KDT.DFS(KDT.root);			//We print after modifications
		*/

		/*KDT.borrar();
		KDT.DFS(KDT.root);			//We print after modifications
		*/

		/*KDT.update();
		KDT.DFS(KDT.root);			//We print after modifications
		*/

	}


	return 0;
}

