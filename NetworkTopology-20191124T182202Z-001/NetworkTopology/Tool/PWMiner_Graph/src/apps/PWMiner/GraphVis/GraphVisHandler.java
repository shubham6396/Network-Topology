package apps.PWMiner.GraphVis;

import apps.PWMiner.GraphModel.Edge;
import apps.PWMiner.GraphModel.Node;
import apps.PWMiner.GraphModel.PWPlainGraph;
import apps.PWMiner.common.Define;
import prefuse.controls.*;
//import apps.PWMiner.GraphModel.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.FileOutputStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Pattern;

/**
 * @author Qiong Cheng
 *
 * Singleton class for handling pathway's image
 */
public class GraphVisHandler extends prefuse.Display {
    //protected static Logger logger = Logger.getLogger(GraphVisHandler.class);

    private static boolean init = false;
    static GraphVisHandler imageLib = new GraphVisHandler();

    private static final String hover = "hover";
    public static final String GRAPH = "graph";
    private static final String nodes = "graph.nodes";

	private static void runProjectCode(String location) throws Exception {

		
		PWPlainGraph testGraph = new PWPlainGraph();
		String currentDirectory = System.getProperty("user.dir");
                String data_location = currentDirectory + "/data/";
		testGraph = testGraph.readGraph(data_location+location+".grp", true);

		int V = testGraph.getNodesNum();
		System.out.println("Number of Vertices  = " + V);
                System.out.println("Vertices->");
		for(int i = 0; i < V; i++){
			System.out.println(testGraph.getNode(i));
		}
                System.out.println(); 
		int E = testGraph.getEdgesNum();
                System.out.println("Number of Edges  = " + E);
                System.out.println("Edges->");        
		for(int i = 0; i < E; i++){
			System.out.println(testGraph.getEdge(i));
		}

                System.out.println();
		int[] degrees = degreeDist(testGraph);
                int []original_degrees=degrees;
                System.out.println("*************Task 0.1( Degree of a node and its distribution)*************");
		System.out.println("Degrees of each node = " + Arrays.toString(degrees));

		int sort_deg[] = degreeDist(testGraph);
		Arrays.sort(sort_deg);
		int K = sort_deg[sort_deg.length-1];
		float[] probies = new float[K];
		int N = countInRange(sort_deg, V, 0, K);
		for(int i = 0; i < K; i++){
			probies[i] = (float) countInRange(sort_deg, V, i+1, i+1) / N;
		}
                System.out.println("Probability of degrees = " + Arrays.toString(probies));
                System.out.println("Average degree = " + findAverage(degrees));
                System.out.println(); 
                System.out.println("*************Task 0.2(Strength of a node and its distribution)*************");
                //task 0.2
                int[] strength=strengthNode(testGraph);
                System.out.println("Strength by node = " + Arrays.toString(strength));
                Arrays.sort(strength);
                
                
                int distinct_count_of_strength=distinctCount(strength);
                
                float [] spread=new float[distinct_count_of_strength];
               
                spread=calculateSpread(strength, distinct_count_of_strength,N);
                System.out.println("Spread by node = " + Arrays.toString(spread));
                
                double averageStrength=findAverage(strength);
                System.out.println("Average Strength of Network = " + averageStrength);
                
                
                
                float C_avg = 0;
		int clusters = 0;
                int[] clusterNodes = localNodes(testGraph);
                float[] C = new float[degrees.length];
                for(int i = 0; i < C.length; i++){
                    if (degrees[i] > 1) {
                        C[i] = (float) clusterNodes[i] / (original_degrees[i] * (original_degrees[i] - 1));
                        C_avg += C[i];
                        clusters += 1;
                    }
                }
                
                C_avg = C_avg/(V);
                System.out.println(); 
                System.out.println("*************Task 0.4(Characteristic path length of a network)*************");
                System.out.println("Shortest path matrix");
                System.out.println("" + print((Shortest_Path(testGraph))));
                double test_CPL=CPL(testGraph);
                System.out.println("Characteristic path length = " + test_CPL);
                System.out.println();
                System.out.println("*************Task 0.5(Clustering coefficient of a network)*************");
                System.out.println("Clustering coefficient by node = " + Arrays.toString(C));
                System.out.println("Average clustering coefficient = " + C_avg);
                
                    
                PWPlainGraph comparison_graph = new PWPlainGraph();
		comparison_graph = testGraph.readGraph(data_location+"Test/bpw2.grp", true);
                
                double comparsion_CPL=CPL(comparison_graph);
                
		int[] comparison_graph_degree = degreeDist(comparison_graph);
                
                float C_avg_comparion = 0;
		int V_comparison = comparison_graph.getNodesNum();
                int[] clusterNodes_comparison = localNodes(comparison_graph);
                float[] C_comparison= new float[V_comparison];
                for(int i = 0; i < C_comparison.length; i++){
                    if (comparison_graph_degree[i] > 1) {
                        C_comparison[i] = (float) clusterNodes_comparison[i] / (comparison_graph_degree[i] * (comparison_graph_degree[i] - 1));
                        C_avg_comparion += C[i];
                        
                    }
                }
                
                C_avg_comparion = C_avg_comparion/(V_comparison);
                System.out.println(); 
                System.out.println("*************Task 0.6(Small World Property of a network)*************");
                System.out.println("Clustering Coefficent Ratio to random network with same number of edges and vertices = " + C_avg/C_avg_comparion);
                System.out.println("Characterstic Path length Ratio with same number of edges and vertices= "+ test_CPL/comparsion_CPL);    
                if((C_avg/C_avg_comparion)/(test_CPL/comparsion_CPL)>1)
                {
                    System.out.println("Graph is of small-world network type");
                }
                else
                {
                    System.out.println("Graph is not small-world network");
                }
                System.out.println();  
                System.out.println("*************Task 0.8(Mixing behavior of nodes)*************");
                double r=mixingBehaviourOfNodes(testGraph);
                System.out.println("R value of Network = " + r);
	        System.out.println("Note: For some in When r = 'NaN', it means One cannot predict the mixing behavior, it is equivalent to r = 0");
                System.out.println(); 
                System.out.println("*************Task 0.10( Closeness Centrality)*************");
             
                closenessCentrality(testGraph);
                System.out.println("*************END OF TASKS*************");
	}

	// an efficient approach to fin N(k) (number of nodes with a given value of degree k)
	// will be to first sort the array and then using a modified binary search
	// function find two indices, one of first element greater than or equal to lower bound
	// of range and the other of the last element less than or equal to upperbound
	// time for running each query will be O(lnV) and for sorting the array once will be O(VlogV).

	// function to find first index >= x
	static int lowerIndex(int arr[], int n, int x)
	{
		int l = 0, h = n - 1;
		while (l <= h)
		{
			int mid = (l + h) / 2;
			if (arr[mid] >= x)
				h = mid - 1;
			else
				l = mid + 1;
		}
		return l;
	}

	// function to find last index <= y
	static int upperIndex(int arr[], int n, int y)
	{
		int l = 0, h = n - 1;
		while (l <= h)
		{
			int mid = (l + h) / 2;
			if (arr[mid] <= y)
				l = mid + 1;
			else
				h = mid - 1;
		}
		return h;
	}

	// function to count elements within given range
	static int countInRange(int arr[], int n, int x, int y)
	{
		// initialize result
		int count = 0;
		count = upperIndex(arr, n, y) -
				lowerIndex(arr, n, x) + 1;
		return count;
	}

	// a basic degree distribution method
	// brute force, runtime O(V^2 E)

	public static int[] degreeDist(PWPlainGraph graph) throws Exception{
		int[] distrib = new int[graph.getNodesNum()];

		for(int i = 0; i < distrib.length; i++){

			// dist[i] = degree of node i
			// loop through edges & count up how many start or end with node i

			Node nodei = graph.getNode(i);
			for(int j = 0; j < distrib.length; j++){

				if (j == i) continue;
				Node nodej = graph.getNode(j);

				// loop through edges and check if either Eij or Eji is in graph

				boolean startend = false;
				boolean endstart = false;

				// check through edges

				for(int e = 0; e < graph.getEdgesNum(); e++){
					Edge dummy = graph.getEdge(e);
					if (dummy.getStartNode().equals(nodei)
							&& dummy.getEndNode().equals(nodej)) {
						startend = true;
						break;
					}
					if (dummy.getStartNode().equals(nodej)
							&& dummy.getEndNode().equals(nodei)) {
						endstart = true;
						break;
					}
				}
				if(startend || endstart) distrib[i]++;
			}
		}
		return distrib;
	}

	public static int findSum(int[] array) {
		int sum = 0;
		for (int value : array) {
			sum += value;
		}
		return sum;
	}

	public static double findAverage(int[] array) {
		int sum = findSum(array);
		return (double) sum / array.length;
	}

        
        //task 0.2
        public static int[] strengthNode(PWPlainGraph graph) throws Exception
        {
            
            int[] strength = new int[graph.getNodesNum()];
        	for(int i = 0; i < strength.length; i++)
                {             
			Node nodei = graph.getNode(i);
			for(int j = 0; j < strength.length; j++)
                        {
                            int possibleInteractions=0;
				if (j == i) continue;
				Node nodej = graph.getNode(j);
				for(int e = 0; e < graph.getEdgesNum(); e++)
                                {
					Edge edge = graph.getEdge(e);
					if (edge.getStartNode().equals(nodei)&& edge.getEndNode().equals(nodej))
                                        {
                                                possibleInteractions++;	
					}
					if (edge.getStartNode().equals(nodej)&& edge.getEndNode().equals(nodei))
                                        {
                                                 possibleInteractions++;		
					}
				}
                                strength[i]=strength[i]+possibleInteractions;
			}
                        
		}
		return strength;
        }
        
        //to count distinct number of strength count
        public static int distinctCount(int [] strength) throws Exception
        {
            int distinct_count_of_strength=1;
            for(int i=1;i<strength.length;i++)
                {
                    if(strength[i]!=strength[i-1])
                    {
                        distinct_count_of_strength++;
                    }              
                }
            return distinct_count_of_strength;
        }
        
        //to calculate Spread P(s)
        public static float[] calculateSpread(int [] strength,int distinct_count_of_strength,int N)
        {
            float [] spread=new float[distinct_count_of_strength];
             int index=0;
                for(int i=0;i<strength.length;i++)
                {
                    
                    int count=1;
                    for(int j=i+1;j<strength.length;j++)
                    {
                        if(strength[i]==strength[j])
                        {
                            count=count+1;
                        }
                        else
                        {
                            
                            i=j-1;
                            break;
                        }
                    }
                    if(index<distinct_count_of_strength)
                    {
                        spread[index]=(float)count/N;
                    index++;
                    }
                    
                    
                }
                return spread;
        }
        
        
        
        //task 0.8
        public static double mixingBehaviourOfNodes(PWPlainGraph graph) throws Exception
        {
              
				
		  double A = 0 ;
		double B = 0;
		double C = 0;
		double D = 0;
		double r = 0;
			
		
	
		try {
			
			
			int V = graph.getNodesNum(); 
			
			int nodeArrayDegrees[] = new int[V];
			int E = graph.getEdgesNum(); 
		
			double edgeCount  = graph.getEdgesNum();
						
			for(int i= 0; i < edgeCount; i++ ) {
				
				Edge edge = graph.getEdge(i);
							
				nodeArrayDegrees[edge.getStartNode().getIndex()]++;
				nodeArrayDegrees[edge.getEndNode().getIndex()]++;
			}
			
			for(int i= 0; i < edgeCount; i++) {
				
				   Edge edge = graph.getEdge(i);			
				   
				   // A Represents summation for j'k'
				   A = A + ((nodeArrayDegrees[edge.getEndNode().getIndex()]) * (nodeArrayDegrees[edge.getStartNode().getIndex()]));  		   
				  
				   // B holds Summation for 0.5 * (j' +  k')
				   B = B + ( (0.5) * ((nodeArrayDegrees[edge.getEndNode().getIndex()])  +  (nodeArrayDegrees[edge.getStartNode().getIndex()]) ) );
			
				   // C represents summation for 0.5 * (j'^2  + k'^2)
				   C = C + 0.5 * ( Math.pow( nodeArrayDegrees[edge.getEndNode().getIndex()] ,  2 )  +  Math.pow(nodeArrayDegrees[edge.getStartNode().getIndex()], 2)  ) ;
			}	
		
		   double inverseEdgeCount =  (1.0 / edgeCount) ;
								
			A = A * inverseEdgeCount;
			
			//System.out.println("A:" + A); 
	        B = B * inverseEdgeCount;
	        
	        B = Math.pow(B, 2);
	        
	        C = C * inverseEdgeCount;	        

	       // Note: For some input grp file if r value is 'NaN', it means One cannot predict the mixing behavior, it is equivalent to r = 0"
            r = (A - B) / (C - B);                
			
		}catch(ArithmeticException e) {
			
			r = 0;
			
		} 
		catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
			
	     
                return r;
		}
               
	
        }
        
        static HashMap<String, String> map;
	static final Pattern tokenizer = Pattern.compile(" ");
	static int[][] nodes2DMatrix;
	static int edgeCount;
        //task 0.10
        public static void closenessCentrality(PWPlainGraph graph) throws Exception
        {
            try
            {
                map = new HashMap<String, String>();
                int V = graph.getNodesNum();
                	edgeCount = graph.getEdgesNum();
			nodes2DMatrix = new int[V][V];
			
			for(int i = 0; i < V ; i++) {
				
				for(int j = 0; j < V; j ++) {
				
					if(i != j) {
					
						nodes2DMatrix[i][j] = Integer.MAX_VALUE;
					
					}
				}
			}
			
			int edgeCount = graph.getEdgesNum();
			
			for(int i = 0; i < edgeCount; i++ ) {
				
				Edge edge = graph.getEdge(i);
				
				int startNode =  edge.getStartNode().getIndex();
				int edgeNode =   edge.getEndNode().getIndex();
				
				if(map.containsKey(String.valueOf(startNode))) {
					
					String endNodes = map.get(String.valueOf(startNode));
					
					endNodes = endNodes + " " + edgeNode ; 
					
					map.put(String.valueOf(startNode), endNodes);
					
				}else{
					
					map.put(String.valueOf(startNode), String.valueOf(edgeNode));
					
				}
			
				// back 
				if(map.containsKey(String.valueOf(edgeNode))) {
					
					String startNodes = map.get(String.valueOf(edgeNode));
					
					startNodes = startNodes + " " + startNode ; 
					
					map.put(String.valueOf(edgeNode), startNodes);
					
				}else{
					
					map.put(String.valueOf(edgeNode), String.valueOf(startNode));
					
				}
					
			}
					
			for(int i = 0 ; i < V ; i++) {
				
				for(int j = 0; j < V ; j++) {
								
				   DFS(i, j, 0, i);
					
				}			
			}
			
		double[] nodes = new double[V]; 	
		
		for(int i = 0 ; i < V ; i++) {
				
			for(int j = 0; j < V ; j++ ) {
					
				nodes[i] = nodes[i] + nodes2DMatrix[i][j];
								
				}	
			
		    nodes[i] = (double)V/(double)nodes[i]; 
                    int value=i+1;
			System.out.println("Closeness centrality for node"+value+"->"+nodes[i]);
            }
        }
            catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        }
        
        //DFS Method
	public static void DFS(int NodeXIndex, int NodeYIndex, int curr_count, int temp_x) {
		 
		if(curr_count > edgeCount + 1 ) {
			
			return;
		}
		
		 if(map.containsKey(String.valueOf(temp_x))) {
			 
			 curr_count = curr_count + 1;
			  String endNodes  =  map.get(String.valueOf(temp_x));
		
			  String[] endNodesArr = null;
			  
			  if(endNodes.contains(" ")) {
	
				  endNodesArr = tokenizer.split(endNodes);
				  
				 for(String tempEndNode: endNodesArr) {
					 
					 if(tempEndNode.equals(String.valueOf(NodeYIndex))){
						 

						 					 
						 if(curr_count < nodes2DMatrix[NodeXIndex][NodeYIndex]) {
							 
							 nodes2DMatrix[NodeXIndex][NodeYIndex] = curr_count;
						 }
						 
					//	break;
					 }else {
						 
						 DFS(NodeXIndex, NodeYIndex, curr_count, Integer.parseInt(tempEndNode));
						 
					 }					 
				 }
			  }else {				  
				  
				  if(endNodes.equals(String.valueOf(NodeYIndex))) {

						 if(curr_count < nodes2DMatrix[NodeXIndex][NodeYIndex]) {
							 
							 nodes2DMatrix[NodeXIndex][NodeYIndex] = curr_count;
						 }
					  
				  }else {
					  
					  DFS(NodeXIndex, NodeYIndex, curr_count, Integer.parseInt(endNodes));
					  
				  }
			  }
			 
		 }
	}
        
        //task 0.5
        public static int[] localNodes(PWPlainGraph graph) throws Exception{

        int N = graph.getNodesNum();
        int[] C = new int[graph.getNodesNum()];
        for(int i = 0; i < N; i++){

            Node node_i = graph.getNode(i);
            for(int j = 0; j < N; j++){

                if (j == i) continue;
                Node node_j = graph.getNode(j);

                // loop through edges and check if either Eij or Eji is in graph

                boolean startend = false;
                boolean endstart = false;

                // check through edges

                for(int e = 0; e < graph.getEdgesNum(); e++){
                    Edge dummy = graph.getEdge(e);
                    if (dummy.getStartNode().equals(node_i)
                            && dummy.getEndNode().equals(node_j)) {
                        startend = true;
                       
                        break;
                    }
                    if (dummy.getStartNode().equals(node_j)
                            && dummy.getEndNode().equals(node_i)) {
                        endstart = true;
                        break;
                    }
                }

                if (startend || endstart) {

                    for(int k = 0; k < N; k++){
                        if (k == j) continue;
                        if (k == i) continue;
                        Node node_k = graph.getNode(k);

                        boolean flag1 = false;
                        boolean flag2 = false;
                        boolean flag3 = false;
                        boolean flag4 = false;

                        for (int e = 0; e < graph.getEdgesNum(); e++) {
                            Edge ed = graph.getEdge(e);
                            if (ed.getStartNode().equals(node_j)
                                    && ed.getEndNode().equals(node_k)) {
                                flag1 = true;
                                break;
                            }
                            if (ed.getStartNode().equals(node_k)
                                    && ed.getEndNode().equals(node_j)) {
                                flag2 = true;
                                break;
                            }
                        }
                        for (int e = 0; e < graph.getEdgesNum(); e++) {
                            Edge ed = graph.getEdge(e);
                            if (ed.getStartNode().equals(node_i)
                                    && ed.getEndNode().equals(node_k)) {
                                flag3 = true;
                                break;
                            }
                            if (ed.getStartNode().equals(node_k)
                                    && ed.getEndNode().equals(node_i)) {
                                flag4 = true;
                                break;
                            }
                        }
                        if(flag1 || flag2) {
                            if (flag3 || flag4) {
                                C[i]++;
                            }
                        }
                    }
                }
            }
        }
        return C;
    }
        
        //0.4
        public static double CPL(PWPlainGraph G) throws Exception{
            int[][] D = Shortest_Path(G);
            int avg = 0;
            int count = 0;
            for (int i = 0; i < D.length; i++){
                for(int j = 0; j < D[i].length; j++){
                    if(D[i][j] == Integer.MAX_VALUE || D[i][j] == 0) continue;
                    avg += D[i][j];
                    count++;
                }
            }
            return (double)avg / (double)count;
        }
        
        public static String print(int[][] D){
            String s = "";
            
            for(int i = 0; i < D.length; i++){
                for(int j = 0; j < D[i].length; j++){
                    if(D[i][j] == Integer.MAX_VALUE) s += "\u221E ";
                    else s += D[i][j] + " ";
                }
                s += "\n";
            }
            
            return s;
        }
        
        public static int[][] Shortest_Path (PWPlainGraph G) throws Exception {
            int N = G.getNodesNum();
            int[][] D = new int[N][N];
            int i,j,k;
            
            //first block: load D with the **adj matrix** for G
            /* D[i][j] = 0 if i = j
             *           1 if edge IJ is in E
             *           infinity otherwise
             */
            for (i = 0; i < N; i++) {
                Node nodei = G.getNode(i);
                for (j = 0; j < N; j++) {
                    if(i == j) {
                        D[i][j] = 0;
                        continue;
                    }
                    Node nodej = G.getNode(j);
                    for(int e = 0; e < G.getEdgesNum(); e++){
                        Edge edge = G.getEdge(e);
                        if(edge.getStartNode().equals(nodei)
                                && edge.getEndNode().equals(nodej)){
                            D[i][j] = 1;
                            break;
                        }
                    }
                    if(D[i][j] == 0) D[i][j] = Integer.MAX_VALUE; 
                    //^infinite distance for nodes not connected
                }
            }
            
            
            //second block: use Floyd's alg to compute the shortest path
            for (k = 0; k < N; k++) {
                for (i = 0; i < N; i++) {
                    for (j = 0; j < N; j++) {
                        int sum = extendSum(D[i][k],  D[k][j]);
                        if (sum < D[i][j]) {
                            D[i][j] = sum;
                        } 
                    } 
                } 
            }
            
            return D;
        }
        
        public static int extendSum(int x, int y){
            //using integer.max_value to simulate infinity
            if(x == Integer.MAX_VALUE || y == Integer.MAX_VALUE) return Integer.MAX_VALUE;
            else return x + y;
        }
        
	private GraphVisHandler() {
		// TODO Auto-generated constructor stub
	}

	public static void init() {
        if (init) {
            return;
        }
        init = true;
    }

	public static PWGraphSearchViz getImage(String basedir, String graphFile, boolean isOrdered, String imageName, int imageL, int imageW){
		PWGraphSearchViz viz=null;
		init();
    	try{
    		viz = new PWGraphSearchViz(basedir  + "/" + Define.DATA_DIR + "/" + graphFile, isOrdered, imageL, imageW, 0);
    	    //ZoomToFitAction fitaction = new ZoomToFitAction(viz.getVisualization());
    	    //fitaction.run(1.0);
    		System.out.println("output " + (viz==null? "T" : "F"));
    		if (viz.isVisible()){
	    	    while ( ! viz.isSaveImageEnable() ){

	    	    }

	    	    System.out.println("output = " + basedir +"/" + Define.GRAPH_VIZ_DIR + "/" + imageName + Define.VIZ_G_SUFFIX);
	    	    FileOutputStream fout = new FileOutputStream(basedir +"/" + Define.GRAPH_VIZ_DIR + "/" + imageName + Define.VIZ_G_SUFFIX);
		    	viz.saveImage(fout, "PNG", 1.0);//1.0/viz.getScale()
		    	fout.close();
    		}else
    			viz=null;
	    }catch (Exception ex){
	    	ex.printStackTrace();
	    }
	    return viz;
    }

	public static PWGraphSearchViz getImageByType(String basedir, String graphFile, boolean isOrdered, String imageName, int imageL, int imageW, int type){
		PWGraphSearchViz viz=null;
		init();
    	try{
    		viz = new PWGraphSearchViz(basedir + "/" + Define.DATA_DIR + "/" + graphFile, isOrdered, imageL, imageW, type);
     		if (! viz.isVisible())
    			viz=null;
	    }catch (Exception ex){
	    	ex.printStackTrace();
	    }
	    return viz;
    }

	public static void saveImageByType(PWGraphSearchViz viz, String basedir, String graphFile1, boolean isOrdered, String imageName, int imageL, int imageW, int type){
		init();
    	try{
    		if (viz.isVisible()){
	    	    while ( ! viz.isSaveImageEnable() ){

	    	    }

	    	    System.out.println("output = " + basedir  +"/" + Define.GRAPH_VIZ_DIR + "/" + imageName + Define.VIZ_G_SUFFIX);
	    	    FileOutputStream fout = new FileOutputStream(basedir +"/" + Define.GRAPH_VIZ_DIR + "/" + imageName + Define.VIZ_G_SUFFIX);
		    	viz.saveImage(fout, "PNG", 1.0);//1.0/viz.getScale()
		    	fout.close();
    		}
	    }catch (Exception ex){
	    	ex.printStackTrace();
	    }
    }

    public static void main(String[] args) throws Exception {
    	//System.setProperty("java.awt.headless","true");
final String pw = "Test/-beta--D-glucuronide_degradation";
		runProjectCode(pw);
    	//ConfigFile.init();
    	//String pw = "Escherichia_coli_K12/superpathway_of_threonine_metabolism";
    	
		final String graphFile = pw + ".grp";
    	final PWGraphSearchViz viz = GraphVisHandler.getImageByType(System.getProperty("user.dir"), graphFile, true, pw, 400, 400, CommonDef.ISMAPPING);
        JFrame frame = new JFrame("Graph Visualization");

        //main display controls
        viz.setSize(425,425);
        viz.pan(5, 5);
        viz.addControlListener(new DragControl());
        viz.addControlListener(new PanControl());
        viz.addControlListener(new ZoomControl());
        viz.addControlListener(new WheelZoomControl());
        viz.addControlListener(new ZoomToFitControl());

        //frame.getContentPane().add(viz);

        // create a new JSplitPane to present the interface
    	JPanel fpanel = new JPanel();
    	JButton saveButton = new JButton(CommonDef.SAVE, CommonDef.createImageIcon("b1.gif", CommonDef.SAVE));
    	saveButton.setPressedIcon(CommonDef.createImageIcon("b1d.gif", CommonDef.SAVE));
    	saveButton.setRolloverIcon(CommonDef.createImageIcon("b1d.gif", CommonDef.SAVE));
    	saveButton.setDisabledIcon(CommonDef.createImageIcon("b1.gif", CommonDef.SAVE));
    	saveButton.setMargin(new Insets(0,0,0,0));
    	saveButton.addActionListener(new ActionListener(){

			public void actionPerformed(ActionEvent arg0) {
				// TODO Auto-generated method stub
				GraphVisHandler.saveImageByType(viz, System.getProperty("user.dir") , graphFile, true, pw, 400, 400, CommonDef.ISMAPPING);
				System.out.println("---------");

				//Upload

			}

    	});
    	fpanel.add(saveButton);

    	JSplitPane split = new JSplitPane(JSplitPane.VERTICAL_SPLIT, viz, fpanel);
    	split.setOneTouchExpandable(true);
    	split.setContinuousLayout(false);
    	split.setDividerLocation(400);

    	// now we run our action list
    	viz.getVisualization().run("draw");

    	frame.getContentPane().add(split);
        frame.setSize(viz.getWidth()+10,viz.getHeight()+160);

        frame.pack();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().setLayout(new BoxLayout(frame.getContentPane(), BoxLayout.Y_AXIS));

        frame.setVisible(true);

    }


}
