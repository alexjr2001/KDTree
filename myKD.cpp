/*
PROYECTO: KD-TREE
AUTOR: ALEXANDER SEBASTIÁN GÓMEZ DEL CARPIO
FECHA: 28/04/2022

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
struct Node {			//Estructura nodo
	Node<T, dim>* right = NULL;			//Punteros a los costados
	Node<T, dim>* left = NULL;
	vector<T> val;
	int nivel = 0;			//Nivel en el que se encuentra del árbol
	Node() {}
	~Node() {}
};

template<class T, int dim>		//Template T del tipo de dato y dim del numero de dimensiones
struct KDTree {				
	Node<T, dim>* root = NULL;
	KDTree() {}
	~KDTree() {}

	void DFS(Node<T, dim>*& ptr) {		//Sirve para imprimir el árbol en cualquier dimensión
		if (ptr == nullptr) {
			return;
		}
		cout << endl << ptr->val[0] << " " << ptr->val[1] << " Nivel: " << ptr->nivel;
		DFS(ptr->right);
		DFS(ptr->left);
	}

	double distance(Node<T,dim>* a,vector<T> b) {		//Para hallar la distancia entre un nodo y una coordenada.
		double dist = 0;
		double suma = 0;
		for (int i = 0; i < dim; i++)
		{
			suma += pow(a->val[i] - b[i], 2);
		}
		dist = sqrt(suma);
		return dist;
	}

	void rangeQuery2(vector<T> query_i,vector<T> query_f,Node<T,dim>* ptr, vector<Node<T,dim>*> &respuesta) {		//Función recursiva para hallar todos los puntos del árbol dentro de un rango.
		if (ptr) {
			int discriminador = ptr->nivel % dim;
			if (query_i[discriminador] <= ptr->val[discriminador] && ptr->val[discriminador] <= query_f[discriminador]) {		//Si la coordenada se encuentra dentro del rango seguimos buscando para la izquierda y derecha
				bool contiene = true;
				for (int i = 0; contiene && i < dim; i++) {							//Comprobar si todas sus dimensiones están dentro del rango
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
			else if (query_i[discriminador] > ptr->val[discriminador]) {		//Si no está en el rango buscamos un punto que si lo esté
				rangeQuery2(query_i, query_f, ptr->right, respuesta);
			}
			else {
				rangeQuery2(query_i, query_f, ptr->left, respuesta);
			}
		}
		return;
	}

	void rangeQuery() {		//Su función es obtener el rango de búsqueda y luego llamar a rangeQuery2() para devolver los puntos en ese rango.
		vector<T> query_i;
		vector<T> query_f;
		vector<Node<T,dim>*> respuesta;
		T v1,v2;
		cout << "\nIngresar intervalos (de menor a mayor) de busqueda\n";
		for (int i = 0; i < dim; i++)
		{
			cout << "\tPara dimension "<< i << ": ";
			cin >> v1 >> v2;
			query_i.push_back(v1); //v1 es el menor
			query_f.push_back(v2); //v2 es mayor
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

	void NN2(vector<T> pto,Node<T,dim>* now,Node<T,dim>*& candidate, int depth, double& bestD=-1) {			//Encuentra el punto más cercano al punto ingresado en NN con recursividad.
		if (!now) {
			return;
		}
		if (bestD==-1 || distance(now, pto) < bestD) {			//Si la distancia actual es mejor de la ya guardada la reemplazamos
			bestD = distance(now, pto);
			candidate = now;
		}

		int axis = depth % dim;
		bool right = false;
		if (pto[axis] < now->val[axis]) {				//Vemos si el punto está al lado derecho o izquierdo
			right = false;
			NN2(pto,now->left,candidate,depth++,bestD);
		}
		else {
			right = true;
			NN2(pto,now->right, candidate, depth++, bestD);
		}
		if (fabs(now->val[axis] - pto[axis]) < candidate->val[axis]) {			//Si la distancia actual haciendo una circunferencia supera el lado a ignorar
			if (right) {
				NN2(pto,now->left, candidate, depth++, bestD);
			}
			else {
				NN2(pto,now->right, candidate, depth++,bestD);
			}
		}
	}

	void NN() {		//Obtiene un punto no existente en el árbol y llama a NN().
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
	void KNN2(vector<T> query, vector<pair<double,Node<T, dim>*>> &ans, Node<T, dim>* ptr) {		// Recursivamente recorre el árbol y halla los K puntos más cercanos al punto ingresado en KNN.
		if (!ptr) {
			return;
		}
		ans.push_back(make_pair(distance(ptr,query),ptr));
		KNN2(query,ans,ptr->right);		
		KNN2(query,ans,ptr->left);
	}

	void KNN() {		// Obtiene un punto no existente en el árbol y llama a KNN().
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

		sort(ans.begin(),ans.end());		//Ordena todas las distancias al punto ingresado
		
		for (int i = 0; i < ans.size() && i < K; i++) {
			cout << "\n\t";
			for (int j = 0; j < dim; j++) {
				cout<<ans[i].second->val[j]<<" ";
			}
		}

		return;
	}

	Node<T, dim>* search(Node<T, dim>* root2, Node<T, dim>*& preptr, vector<T>& query, int n) {			//Busca un nodo en el árbol
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
			if (ptr->val[nivel % dim] > query[nivel % dim]) {   // A la izquierda si es menor el valor a buscar
				preptr = ptr;
				ptr = ptr->left;
			}
			else if (ptr->val[nivel % dim] < query[nivel % dim]) {   //A la derecha si es mayor
				preptr = ptr;
				ptr = ptr->right;
			}
			else {						//Si son iguales 
				resultado = ptr;
				for (int i = 0; resultado && i < dim; i++) {		//Vemos si es el punto que se está buscando
					if (ptr->val[i] != query[i]) {
						resultado = nullptr;
					}
				}
				if (resultado) return resultado;
				else {												//Si no lo es vamos hacia los dos lados
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

	bool insert() {				//Inserta un nodo en el árbol
		vector<T> query;
		cout << "\nInserta las nuevas coordenadas: ";
		T v;
		for (int i = 0; i < dim; i++) {
			cin >> v;
			query.push_back(v);
		}
		Node<T, dim>* preptr = root;						//Apuntará al nodo un nivel arriba del lugar de inserción del nuevo nodo
		Node<T, dim>* ptr = search(root, preptr, query, 0);			//Si el nodo ya está en el árbol ya no lo insertamos

		Node<T, dim>* nuevo = new Node<T, dim>;
		if (ptr) {
			cout << "El valor ya existe, ingrese otro nodo\n";
			return false;
		}
		else {
			nuevo->val = query;
			if (!preptr) {				//Si el árbol está vacío
				root = nuevo;
			}
			else {
				nuevo->nivel = preptr->nivel + 1;
				if (query[preptr->nivel % dim] > preptr->val[preptr->nivel % dim]) {		//Vemos si lo insertamos a la derecha o a la izquierda
					preptr->right = nuevo;
				}
				else {
					preptr->left = nuevo;
				}
			}
			return true;
		}
	}



	void findMin(Node<T, dim>*& min, Node<T, dim>* tmp, int d, Node<T, dim>*& preptr, Node<T, dim>* pretmp) {			//Encuentra el mínimo valor en una dimension específica en un subárbol
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

	void borrar2(Node<T, dim>* ptr, Node<T, dim>* preptr, Node<T, dim>* min = nullptr) {	//Función recursiva para borrar un nodo
		if (ptr) {	//Es que se encontró
			if (ptr->right) {				//Si no es hoja y tiene hijo a la derecha
				min = ptr->right;
				preptr = ptr;
				findMin(min, ptr->right, ptr->nivel % dim, preptr, ptr);		//Debemos encontra el mínimo en la dimensión discriminante y en su subárbol
				ptr->val = min->val;
				borrar2(min, preptr);
			}
			else if (ptr->left) {			//Tiene hijo a la izquierda que luego serán los hijos a la derecha
				min = ptr->left;
				preptr = ptr;
				findMin(min, ptr->left, ptr->nivel % dim, preptr, ptr);
				ptr->val = min->val;
				ptr->right = ptr->left;
				ptr->left = nullptr;
				borrar2(min, preptr);
			}
			else {							//Si es hoja solo se elimina
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
	bool borrar() {				//Obtiene el valor a borrar
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

	void update() {			//Actualiza un valor en el árbol
		vector<T> nuevas;
		vector<T> query;
		if (borrar()) {
			while(!insert()){}
		}
		return;
	}
};


template<class T, int dim>
struct compare						//Ayuda a la función sort de un vector a comparar para poder ordenarlo
{
	int d = 0; //Dimension a comparar
	compare(int _d) :d(_d) {}
	inline bool operator() (const Node<T, dim>* struct1, const Node<T, dim>* struct2)
	{
		return (struct1->val[d] < struct2->val[d]);
	}
};


template<class T, int dim>
void sortV(vector<Node<T, dim>*>& v, int d) {   //Ordena en una dimension
	std::sort(v.begin(), v.end(), compare<T, dim>(d));
	return;
}

template<class T, int dim>
void armarTree(KDTree<T, dim>& KD, vector<Node<T, dim>*>& v, Node<T, dim>*& ptr, int d) {		//Con datos ingresados al inicio se arma el árbol de inicio
	int size = v.size();
	int i_nodo = size / 2;
	if (size == 0) {
		return;
	}
	sortV(v, d % dim);		//Ordena el vector de nodos en la dimensión correspondiente
	ptr = v[i_nodo];
	v[i_nodo]->nivel = d;
	vector<Node<T, dim>*> v_left;
	copy(v.begin(), v.begin() + i_nodo, back_inserter(v_left));		//Cortamos el vector en dos para llamar armarTree() para su respectiva construcción
	vector<Node<T, dim>*> v_right;
	copy(v.begin() + i_nodo + 1, v.end(), back_inserter(v_right));

	armarTree(KD, v_left, ptr->left, d + 1);
	armarTree(KD, v_right, ptr->right, d + 1);

	return;
}


int main() {
	const int n_nodos = 8;	//Nro de nodos
	const int d = 2;	//dimension
	vector<Node<int, d>*> v_inicio;			//Almacena los nodos del inicio
	Node<int, d>* nodos = new Node<int, d>[n_nodos]; //Creando nodos de inicio
	KDTree<int, d> KDT;

	//INSTANCIA
	int instancia[n_nodos][d] = {{7,7},{1,2},{3,4},{6,9},{9,7},{18,10},{7,2},{1,10}};
	for (int i = 0; i < n_nodos; i++) {
		v_inicio.push_back(&nodos[i]);
		for (int j = 0; j < d; j++) {
			v_inicio[i]->val.push_back(instancia[i][j]);
		}
	}

	//BUCLES PARA PEDIR VALORES AL USUARIO EN EL TERMINAL//////////////////////////////////////////////
	/*int v;
	for (int i = 0; i < n_nodos; i++) {
		v_inicio.push_back(&nodos[i]);	   //Rellenando vector con coordenadas iniciales
		for (int j = 0; j < d; j++) {		
			cin >> v;
			v_inicio[i]->val.push_back(v);
		}
	}*/
	//////////////////////////////////////////////////////////////////////////////////////

	armarTree(KDT, v_inicio, KDT.root, 0);		//Armamos el árbol

	KDT.DFS(KDT.root);			//Imprimimos el árbol inicial

	int tests = 1;			//Números de tests
	for(int i=0;i<tests ;i++){	//DESCOMENTAR PARA PROBAR CADA FUNCION----------------

		KDT.KNN();				//Los vecinos más cercanos	

		//KDT.NN();						//El vecino más cercano
		
		
		//KDT.rangeQuery();				 //Los puntos dentro de un rango
		

		/*KDT.insert();					//Borrar e insertar tienen la función buscar dentro
		KDT.DFS(KDT.root);			//Imprimimos luego de modificaciones
		*/

		/*KDT.borrar();
		KDT.DFS(KDT.root);			//Imprimimos luego de modificaciones
		*/

		/*KDT.update();
		KDT.DFS(KDT.root);			//Imprimimos luego de modificaciones
		*/

	}


	return 0;
}

