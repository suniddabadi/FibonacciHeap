import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.Scanner;
import java.util.Stack;

/*This class is does for all Fibonacci heap operations
 * mainly Inserting elelment in the fibonacci heap and remove minimum element 
 */
class FibonacciHeap
{
	private FibonacciHeapNode minNode;
	private int nofnodes;

	// constructor
	public FibonacciHeap()
	{

	}

	public boolean isEmpty() // Returns null if heap is empty
	{
		return minNode == null;
	}

	public void clear() // Resets the heap
	{
		minNode = null;
		nofnodes = 0;
	}

	// insert functions inserts the node in the heap in O(1) time.
	public FibonacciHeapNode insert(FibonacciHeapNode node, double key)
	{
		node.key = key;
		if (minNode != null)
		{
			node.left = minNode;
			node.right = minNode.right;
			minNode.right = node;
			node.right.left = node;

			if (key < minNode.key)
			{
				minNode = node;
			}

		}
		else
		{
			minNode = node;
		}

		nofnodes++;
		return node;
	}

	public FibonacciHeapNode min() // This function returns Minimum node in the
									// heap.
	{
		return minNode;
	}

	public int size() // This function returns the current size of the heap
	{
		return nofnodes;
	}

	/*
	 * This function is a helper function for consolidate function for linking
	 * of node y and node x in the heap
	 */
	protected void link(FibonacciHeapNode y, FibonacciHeapNode x)
	{
		// remove y from root list of heap
		y.left.right = y.right;
		y.right.left = y.left;

		// make y a child of x
		y.parent = x;

		if (x.child == null)
		{
			x.child = y;
			y.right = y;
			y.left = y;
		} else
		{
			y.left = x.child;
			y.right = x.child.right;
			x.child.right = y;
			y.right.left = y;
		}

		// increase degree[x]
		x.degree++;

		// set mark[y] false
		y.mark = false;
	}

	//Declaration of fibonacci heap property
	private static final double oneOverLogPhi = 1.0 / Math.log((1.0 + Math
			.sqrt(5.0)) / 2.0);

	public FibonacciHeapNode removeMin()
	{
		FibonacciHeapNode z = minNode;

		if (z != null)
		{
			int numKids = z.degree;
			FibonacciHeapNode x = z.child;
			FibonacciHeapNode tempRight;

			// for each child of z do...
			while (numKids > 0)
			{
				tempRight = x.right;

				// remove x from child list
				x.left.right = x.right;
				x.right.left = x.left;

				// add x to root list of heap
				x.left = minNode;
				x.right = minNode.right;
				minNode.right = x;
				x.right.left = x;

				// set parent[x] to null
				x.parent = null;
				x = tempRight;
				numKids--;
			}

			// remove z from root list of heap
			z.left.right = z.right;
			z.right.left = z.left;

			if (z == z.right)
			{
				minNode = null;
			} else
			{
				minNode = z.right;
				combine();
			}

			// decrement size of heap
			nofnodes--;
		}

		return z;
	}

	protected void combine()
	{
		int arraySize = ((int) Math.floor(Math.log(nofnodes) * oneOverLogPhi)) + 1;

		ArrayList<FibonacciHeapNode> array = new ArrayList<FibonacciHeapNode>(
				arraySize);

		// Initialize degree array
		for (int i = 0; i < arraySize; i++)
		{
			array.add(null);
		}

		// Find the number of root nodes.
		int numRoots = 0;
		FibonacciHeapNode x = minNode;

		if (x != null)
		{
			numRoots++;
			x = x.right;

			while (x != minNode)
			{
				numRoots++;
				x = x.right;
			}
		}

		// For each node in root list do...
		while (numRoots > 0)
		{
			// Access this node's degree..
			int d = x.degree;
			FibonacciHeapNode next = x.right;

			// ..and see if there's another of the same degree.
			for (;;)
			{
				FibonacciHeapNode y = array.get(d);
				if (y == null)
				{
					// Nope.
					break;
				}

				// There is, make one of the nodes a child of the other.
				// Do this based on the key value.
				if (x.key > y.key)
				{
					FibonacciHeapNode temp = y;
					y = x;
					x = temp;
				}

				// FibonacciHeapNode<T> y disappears from root list.
				link(y, x);

				// We've handled this degree, go to next one.
				array.set(d, null);
				d++;
			}

			// Save this node for later when we might encounter another
			// of the same degree.
			array.set(d, x);

			// Move forward through list.
			x = next;
			numRoots--;
		}

		// Set min to null (effectively losing the root list) and
		// reconstruct the root list from the array entries in array[].
		minNode = null;

		for (int i = 0; i < arraySize; i++)
		{
			FibonacciHeapNode y = array.get(i);
			if (y == null)
			{
				continue;
			}

			// We've got a live one, add it to root list.
			if (minNode != null)
			{
				// First remove node from root list.
				y.left.right = y.right;
				y.right.left = y.left;

				// Now add to root list, again.
				y.left = minNode;
				y.right = minNode.right;
				minNode.right = y;
				y.right.left = y;

				// Check if this is a new min.
				if (y.key < minNode.key)
				{
					minNode = y;
				}
			} else
			{
				minNode = y;
			}
		}
	}

	protected void cascadeCut(FibonacciHeapNode y)
	{
		FibonacciHeapNode z = y.parent;

		// if parent exists...
		if (z != null)
		{
			// if y is unmarked, set it marked
			if (!y.mark)
			{
				y.mark = true;
			} else
			{
				// it's marked, cut it from parent
				cut(y, z);

				// cut its parent as well
				cascadeCut(z);
			}
		}
	}

	public void decreaseKey(FibonacciHeapNode x, double k)
	{
		if (k > x.key)
		{
			return;
		}

		x.key = k;

		FibonacciHeapNode y = x.parent;

		if ((y != null) && (x.key < y.key))
		{
			cut(x, y);
			cascadeCut(y);
		}

		if (x.key < minNode.key)
		{
			minNode = x;
		}
	}

	protected void cut(FibonacciHeapNode x, FibonacciHeapNode y)
	{
		// remove x from childlist of y and decrement degree[y]
		x.left.right = x.right;
		x.right.left = x.left;
		y.degree--;

		// reset y.child if necessary
		if (y.child == x)
		{
			y.child = x.right;
		}

		if (y.degree == 0)
		{
			y.child = null;
		}
   
		// add x to root list of heap
		x.left = minNode;
		x.right = minNode.right;
		minNode.right = x;
		x.right.left = x;

		// set parent[x] to nil
		x.parent = null;

		// set mark[x] to false
		x.mark = false;
	}

	//Fibonacci Heap implementation for the randomly generated graph 
	public void Dijk_Fh(Random_Graph graph, int source)
	{
		// Decalaration of the local variables
		Long start = System.currentTimeMillis();
		int[] cost = new int[graph.noofvertices];
		for (int i = 1; i < cost.length; i++)
		{
			cost[i] = 0;
		}
		FibonacciHeapNode[] refofHeap = new FibonacciHeapNode[graph.noofvertices];
		ArrayList<Integer> visited = new ArrayList<>();
		FibonacciHeapNode node = new FibonacciHeapNode(source, 0);
		refofHeap[source] = insert(node, 0);
		// Inserting 1st node in the heap
		for (int i = 0; i < graph.noofvertices; i++)
		{

			node = new FibonacciHeapNode(i, 99999);
			// Inserting remaining nodes in on the heap
			refofHeap[i] = insert(node, 99999);
		}

		FibonacciHeapNode x = removeMin();
		// Removing minimum element from a tree
		cost[source] = 0;
		visited.add(x.value);
		for (Edge j = graph.adjecencylist[source].adjList; j != null; j = j.next)
		{
			if (!visited.contains(j.vertex_number))
			{

				decreaseKey(refofHeap[j.vertex_number], j.cost);
			}
		}

		for (int i = 1; i < graph.noofvertices; i++)
		{
			x = removeMin();
			cost[x.value] = cost[x.value] + (int) x.key;

			visited.add(x.value);
			for (Edge j = graph.adjecencylist[x.value].adjList; j != null; j = j.next)
			{
				if (!visited.contains(j.vertex_number))
				{
					decreaseKey(refofHeap[j.vertex_number], j.cost
							+ cost[x.value]);
				}
			}
		}

		long stop = System.currentTimeMillis();

		System.out.println("Time required for an Execution:   "
				+ (stop - start) + " Miliseconds"); // Printing time required
	}

	//Fibonacci Heap implementation for the graph being read from file
	public void Dijk_Fh(Dijkstra graph, int source)
	{
		// Decalaration of the local variables
		Long start = System.currentTimeMillis();
		int[] cost = new int[graph.noofvertices];
		for (int i = 1; i < cost.length; i++)
		{
			cost[i] = 0;
		}
		FibonacciHeapNode[] refofHeap = new FibonacciHeapNode[graph.noofvertices];
		ArrayList<Integer> visited = new ArrayList<>();
		FibonacciHeapNode node = new FibonacciHeapNode(source, 0);
		refofHeap[source] = insert(node, 0); // Inserting 1st node in the heap
		for (int i = 0; i < graph.noofvertices; i++)
		{

			node = new FibonacciHeapNode(i, 99999);// Inserting remaining nodes
													// in on the heap
			refofHeap[i] = insert(node, 99999);
		}

		FibonacciHeapNode x = removeMin();// Removing minimum element from a
											// tree
		cost[source] = 0;
		visited.add(x.value);
		for (Edge j = graph.adjecencylist[source].adjList; j != null; j = j.next)
		{
			if (!visited.contains(j.vertex_number))
			{

				decreaseKey(refofHeap[j.vertex_number], j.cost);

			}
		}

		for (int i = 1; i < graph.noofvertices; i++)
		{
			x = removeMin();
			cost[x.value] = cost[x.value] + (int) x.key;
			// System.out.println(""+ x.data+"  "+vertex_Array[x.data]);
			visited.add(x.value);
			for (Edge j = graph.adjecencylist[x.value].adjList; j != null; j = j.next)
			{
				if (!visited.contains(j.vertex_number))
				{
					decreaseKey(refofHeap[j.vertex_number], j.cost
							+ cost[x.value]);
				}
			}
		}

		long stop = System.currentTimeMillis();
		for (int i = 0; i < cost.length; i++)
		{
			System.out.println(cost[i]); // Printing final cost
		}
		System.out.println("Time required for an Execution:   "
				+ (stop - start) + " Miliseconds"); // Printing time required
	}
}

/*
 * Implements a node of the Fibonacci heap. It holds the information necessary
 * for maintaining the structure of the heap. It also holds the reference to the
 * key value.
 */
class FibonacciHeapNode
{

	int value;

	FibonacciHeapNode child;
	FibonacciHeapNode left;
	FibonacciHeapNode right;
	FibonacciHeapNode parent;
	boolean mark;
	double key;
	int degree;

	public FibonacciHeapNode(int data, double key)
	{

		right = this;
		left = this;
		this.value = data;
		this.key = key;
	}

	public final double getKey()
	{
		return key;
	}

	public final int getData()
	{
		return value;
	}

}

/* Implements the vertex of the graph */
class Vertex
{

	int number;
	Edge adjList;
	boolean visited;

	public Vertex(int number, Edge adjList, boolean visited)
	{

		this.number = number;
		this.adjList = adjList;
		this.visited = visited;
	}

	public Vertex()
	{

		this.number = 0;
		this.adjList = null;
		this.visited = false;

	}

}

/* Implements Edges of the graph */
class Edge
{

	int vertex_number;
	int cost;
	Edge next;

	public Edge(int vertex_number, int cost, Edge next)
	{

		this.vertex_number = vertex_number;
		this.cost = cost;
		this.next = next;
	}

	public Edge()
	{
		this.vertex_number = 0;
		this.cost = 0;
		this.next = null;
	}
}

/*
 * This class creates the random graph
 */
class Random_Graph
{
	int noofvertices;
	int noofedges;
	Vertex[] adjecencylist;
	int density;

	public Random_Graph(int vertices, double d)
	{
		int max_edges = (vertices * (vertices - 1)) / 2;
		float divide_density = (float) (d / 100f);
		double possible = Math.ceil((max_edges) * divide_density);

		noofvertices = vertices;
		noofedges = (int) possible;
	}
	/*
	 * This function creates a random graph.and returns an instance of it
	 */
	public void createGraph()
	{
		Random generator = new Random();
		Random cost_generator = new Random();
		int vertex1, vertex2, edge_cost;
		int edge_count = 0;
		boolean b[][] = new boolean[noofvertices][noofvertices];

		adjecencylist = new Vertex[noofvertices];
		for (int i = 0; i < noofvertices; i++)
		{
			adjecencylist[i] = new Vertex(i, null, false);
		}
		vertex1 = generator.nextInt(noofvertices);
		vertex2 = generator.nextInt(noofvertices);
		while (vertex1 == vertex2)
		{
			vertex2 = generator.nextInt(noofvertices);
		}
		edge_cost = cost_generator.nextInt(1000) + 1;
		if (!b[vertex1][vertex2] && !b[vertex2][vertex1])
		{
			adjecencylist[vertex1].adjList = new Edge(vertex2, edge_cost,
					adjecencylist[vertex1].adjList);
			adjecencylist[vertex2].adjList = new Edge(vertex1, edge_cost,
					adjecencylist[vertex2].adjList);
			b[vertex1][vertex2] = true;
			b[vertex2][vertex1] = true;
			edge_count++;
		}

		while (edge_count < noofedges || !depthFisrtSearch(this))
		{
			vertex1 = generator.nextInt(noofvertices);
			vertex2 = generator.nextInt(noofvertices);
			while (vertex1 == vertex2)
			{

				vertex2 = generator.nextInt(noofvertices);
			}
			edge_cost = cost_generator.nextInt(1000) + 1;
			if (!b[vertex1][vertex2] && !b[vertex2][vertex1])
			{
				adjecencylist[vertex1].adjList = new Edge(vertex2, edge_cost,
						adjecencylist[vertex1].adjList);
				adjecencylist[vertex2].adjList = new Edge(vertex1, edge_cost,
						adjecencylist[vertex2].adjList);
				b[vertex1][vertex2] = true;
				b[vertex2][vertex1] = true;
				edge_count++;
			}

		}
		System.out.println("Total Edges are:" + edge_count);

	}

	/*
	 * returns true if the graph is connected
	 */

	public boolean depthFisrtSearch(Dijkstra Graph_for_dfs)
	{
		Stack<Integer> stack = new Stack<Integer>();
		int[] visitedVertices = new int[Graph_for_dfs.noofvertices];
		for (int i = 0; i < visitedVertices.length; i++)
		{
			visitedVertices[i] = -1;
		}

		stack.push(0);
		while (!stack.isEmpty())
		{
			int u = stack.pop();
			if (visitedVertices[u] == -1)
			{
				visitedVertices[u] = 1;
			}
			for (Edge iterator = Graph_for_dfs.adjecencylist[u].adjList; iterator != null; iterator = iterator.next)
			{
				if (visitedVertices[iterator.vertex_number] == -1)
				{
					stack.push(iterator.vertex_number);
				}
			}

		}
		int count = 0;
		for (int i = 0; i < visitedVertices.length; i++)
		{
			if (visitedVertices[i] == 1)
			{
				count++;
			}
		}
		if (count == Graph_for_dfs.noofvertices)
		{
			// System.out.println("Graph is connected!!");
			return true;
		} else
		{
			System.out.println("Graph is not conneceted by "
					+ (Graph_for_dfs.noofvertices - count) + " Nodes");
			return false;
		}
	}

	/*
	 * /* This function returns true if the graph is connected Difference is
	 * this dfs works for random generated graph
	 */
	public boolean depthFisrtSearch(Random_Graph Graph_for_dfs)
	{
		Stack<Integer> stack = new Stack<Integer>();
		int[] visitedvertices = new int[Graph_for_dfs.noofvertices];
		for (int i = 0; i < visitedvertices.length; i++)
		{
			visitedvertices[i] = -1;
		}

		stack.push(0);
		while (!stack.isEmpty())
		{
			int u = stack.pop();
			if (visitedvertices[u] == -1)
			{
				visitedvertices[u] = 1;
			}
			for (Edge iterator = Graph_for_dfs.adjecencylist[u].adjList; iterator != null; iterator = iterator.next)
			{
				if (visitedvertices[iterator.vertex_number] == -1)
				{
					stack.push(iterator.vertex_number);
				}
			}
		}
		
		int count = 0;
		for (int i = 0; i < visitedvertices.length; i++)
		{
			if (visitedvertices[i] == 1)
			{
				count++;
			}
		}
		if (count == Graph_for_dfs.noofvertices)
		{

			return true;
		} else
		{
			return false;
		}
	}

	/* prints the graph */
	public void print()
	{
		for (int i = 0; i < adjecencylist.length; i++)
		{
			for (Edge next = adjecencylist[i].adjList; next != null; next = next.next)
			{
				System.out.println(adjecencylist[i].number + "  "
						+ next.vertex_number + "  " + next.cost);
			}
		}
	}

}

/*
 * This is the main class which reads the input and handles the flow of the
 * program
 */
public class Dijkstra
{

	int noofvertices;
	int noofedges;
	int source;
	Vertex[] adjecencylist;

	public Dijkstra(int no_of_vertices, int no_of_edges, Vertex[] adjecency_list)
	{

		this.noofvertices = no_of_vertices;
		this.noofedges = no_of_edges;
		this.adjecencylist = adjecency_list;
	}

	/* Reads the input file provided and creates the graph */
	public void createGraph(String filename)
	{
		try
		{
			Scanner Reader = null;
			File file1 = new File("C:\\Users\\SUNID\\workspace\\FibonacciHeap\\src\\file.txt");
			Reader = new Scanner(file1);
			this.source = Reader.nextInt();
			this.noofvertices = Reader.nextInt();
			this.noofedges = Reader.nextInt();
			adjecencylist = new Vertex[noofvertices];
			boolean b[][] = new boolean[noofvertices][noofvertices];

			// Reading vertices numbers
			for (int i = 0; i < noofvertices; i++)
			{
				adjecencylist[i] = new Vertex(i, null, false);
			}
			while (Reader.hasNext())
			{
				int vertex1 = Reader.nextInt();
				int vertex2 = Reader.nextInt();
				int edge_cost = Reader.nextInt();
				if (!b[vertex1][vertex2] && !b[vertex2][vertex1])
				{
					adjecencylist[vertex1].adjList = new Edge(vertex2,
							edge_cost, adjecencylist[vertex1].adjList);
					adjecencylist[vertex2].adjList = new Edge(vertex1,
							edge_cost, adjecencylist[vertex2].adjList);
					b[vertex1][vertex2] = true;
					b[vertex2][vertex1] = true;
				}
			}
			Reader.close();

		} catch (Exception e)
		{
			// TODO: handle exception
			System.out.println("Exception occured!!!");
		}

	}

	public Dijkstra()
	{
		this.noofvertices = 0;
		this.noofedges = 0;
		this.adjecencylist = null;

	}

	/* Prints the graph */
	public void print()
	{
		for (int i = 0; i < adjecencylist.length; i++)
		{
			for (Edge next = adjecencylist[i].adjList; next != null; next = next.next)
			{
				System.out.println(adjecencylist[i].number + "  "
						+ next.vertex_number + "  " + next.cost);
			}
		}
	}

	public static void main(String[] args)
	{

		if (args.length < 0)
		{
			System.out.println("Please Provide the input!!");
		} else
		{
			if (args.length == 2)
			{
				String filename = args[1];
				if (args[0].equals("-s"))
				{
					System.out.println("Output for " + args[1]
							+ " as Simple scheme");
					System.out.println("");
					Dijkstra graph = new Dijkstra();
					graph.createGraph(filename);

					new Userfile(graph, graph.source);

				} else
				{
					System.out.println("Output for " + args[1]
							+ "Fibonacci scheme");
					System.out.println("");

					Dijkstra graph = new Dijkstra();
					graph.createGraph(filename);
					FibonacciHeap fh = new FibonacciHeap();
					fh.Dijk_Fh(graph, graph.source);
				}
			} else if (args.length == 4)
			{
				int n = Integer.parseInt(args[1]);
				double d = Double.parseDouble(args[2]);
				int x = Integer.parseInt(args[3]);
				Random_Graph rg = new Random_Graph(n, d);
				rg.createGraph();
				new Userfile(rg, x);
				FibonacciHeap fh = new FibonacciHeap();
				fh.Dijk_Fh(rg, x);

			}
		}
	}
}

class Userfile
{
	Dijkstra graph;
	Random_Graph Rgraph;
	int i = 0, j = 0, k = 0;

	int min = 99999999;
	int[] min_cost;
	int counter = 0;
	int kkkk = 0;
	int final_edges[][];
	ArrayList<Integer> visited;
	int[] distance;
	List<Integer> Settled_Nodes = new ArrayList<Integer>();
	List<Integer> UnSettled_Nodes = new ArrayList<Integer>();

	// calculate Dijkstra for random graph
	public Userfile(Random_Graph graph, int source)
	{

		long start = System.currentTimeMillis();
		this.Rgraph = graph;
		distance = new int[Rgraph.noofvertices];
		for (int i = 0; i < Rgraph.noofvertices; i++)
			distance[i] = min;
		distance[source] = 0;
		UnSettled_Nodes.add(source);
		while (UnSettled_Nodes.size() > 0)
		{

			int v1 = -1;
			for (int i : UnSettled_Nodes)
			{
				if (v1 == -1)
					v1 = i;
				else if (distance[i] < distance[v1])
					v1 = i;
			}

			Settled_Nodes.add(v1);
			UnSettled_Nodes.remove((Integer) v1);

			for (Edge iterator = Rgraph.adjecencylist[v1].adjList; iterator != null; iterator = iterator.next)
			{
				if (Settled_Nodes.contains(iterator.vertex_number))
					continue;
				else
				{
					if (distance[iterator.vertex_number] > (distance[v1] + iterator.cost))
					{
						distance[iterator.vertex_number] = distance[v1]
								+ iterator.cost;
						UnSettled_Nodes.add(iterator.vertex_number);
					}
				}
			}
		}

		long stop = System.currentTimeMillis();
		System.out.println("Simple scheme total time: " + (stop - start)
				+ " Miliseconds");

	}

	//Calculate Dijkstra's for the graph from the file  
	public Userfile(Dijkstra graph, int source)
	{

		long start = System.currentTimeMillis();
		Dijkstra graph1 = graph;
		distance = new int[graph1.noofvertices];
		for (int i = 0; i < graph1.noofvertices; i++)
			distance[i] = min;
		distance[source] = 0;
		UnSettled_Nodes.add(source);
		while (UnSettled_Nodes.size() > 0)
		{
            //checking for the next minimum distance between the unsettled nodes
			int v1 = -1;			
			for (int i : UnSettled_Nodes)
			{
				if (v1 == -1)
					v1 = i;
				else if (distance[i] < distance[v1])
					v1 = i;
			}

			Settled_Nodes.add(v1);
			UnSettled_Nodes.remove((Integer) v1);

			for (Edge iterator = graph1.adjecencylist[v1].adjList; iterator != null; iterator = iterator.next)
			{
				if (Settled_Nodes.contains(iterator.vertex_number))
					continue;
				else
				{
					if (distance[iterator.vertex_number] > (distance[v1] + iterator.cost))
					{
						distance[iterator.vertex_number] = distance[v1]
								+ iterator.cost;
						UnSettled_Nodes.add(iterator.vertex_number);
					}
				}
			}
		}

		long stop = System.currentTimeMillis();
		System.out.println("Simple scheme total time: " + (stop - start)
				+ " Miliseconds");
		System.out.println("The cost from the source to each of the vertices are :");
		for (int i = 0; i < distance.length; i++)
		{
			System.out.println("Cost from source to vertex" + i + " is: " + distance[i]);
		}
		System.out.println("The original graph with the cost was :");
		graph1.print();
	}
}
