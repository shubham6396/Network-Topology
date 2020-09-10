package apps.PWMiner.GraphModel;
/**
 * @author ITCS 6114
 * Pathway Plaintext Graph Handler
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;
import java.util.Vector;

import apps.PWMiner.common.Define;

public class PWPlainGraph{
    //protected static Logger logger = Logger.getLogger(PWPlainGraph.class);

    private Vector nodes ;
    private Vector edges;

    private int type;
    
    public PWPlainGraph(){
    	nodes = new Vector();
    	edges = new Vector();
    }
    
	public void loadGraph(Node[] nodesArr, Edge[] edgesArr) {
		if (edges.size()>0) edges.removeAllElements();
		if (nodes.size()>0) nodes.removeAllElements();
        nodes = new Vector();
        if(nodesArr != null)
        {
            for(int i = 0; i < nodesArr.length; i++)
                nodes.addElement(nodesArr[i]);

        }
        edges = new Vector();
        if(edgesArr != null)
        {
            for(int i = 0; i < edgesArr.length; i++)
                edges.addElement(edgesArr[i]);

        }
	}

	public static PWPlainGraph readGraph(String filename, boolean isOrdered) throws IOException {
        File file = Define.getAbsFile(filename);
        BufferedReader reader = new BufferedReader(new FileReader(file));
        return readGraphByBuffer(reader, isOrdered);
	}        
    
	public static PWPlainGraph readGraphByInputStream(InputStream in, boolean isOrdered) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        return readGraphByBuffer(reader, isOrdered);
	}  
	
	private static PWPlainGraph readGraphByBuffer(BufferedReader reader, boolean isOrdered) throws IOException {
        double x,y;
        int nIndex=0;

        PWPlainGraph graph = new PWPlainGraph();
        boolean hasCord = false;
        for(String line = reader.readLine(); line != null; line = reader.readLine())
        {
            if(line.startsWith("edge"))
                break;
            if(line.indexOf('x') == -1)
                hasCord = false;
            else
                hasCord = true;
            StringTokenizer st = new StringTokenizer(line, " ");
            st.nextToken();
            String label = st.nextToken();

            if(!hasCord)
            {
                x = -1;//(int)(Math.random() * 500D) + 20;
                y = -1;//(int)(Math.random() * 500D) + 20;
            } else
            {
                String cord = st.nextToken();
                x = Double.parseDouble(cord.substring(2));
                cord = st.nextToken();
                y = Integer.parseInt(cord.substring(2));
            }
            Node node = new Node(label, nIndex);
            if (x!=-1 && y!=-1) {
            	node.setX(x);
            	node.setY(y);
            }
            graph.nodes.add(node);
            nIndex++;
        }

        for(String line = reader.readLine(); line != null; line = reader.readLine())
        {
            if(line.startsWith("*"))
                break;
            StringTokenizer st = new StringTokenizer(line, "-");
            int n1 = Integer.parseInt(st.nextToken()) - 1;
            int n2 = Integer.parseInt(st.nextToken()) - 1;
            Edge edge = new Edge((Node)(graph.nodes.get(n1)), (Node)(graph.nodes.get(n2)), isOrdered);
            graph.edges.add(edge);
        }
        
        //System.out.println(" Graph " +  graph.getNodesNum() + " --- " + graph.getEdgesNum());
        reader.close();
        return graph;
	}

	public static PWPlainGraph downloadGraph(String filename, boolean isOrdered) throws IOException {
        double x,y;
        int nIndex=0;

        File file = Define.getAbsFile(filename);
		BufferedReader reader = new BufferedReader(new FileReader(file));
        PWPlainGraph graph = new PWPlainGraph();
        boolean hasCord = false;
        for(String line = reader.readLine(); line != null; line = reader.readLine())
        {
            if(line.startsWith("edge"))
                break;
            if(line.indexOf('x') == -1)
                hasCord = false;
            else
                hasCord = true;
            StringTokenizer st = new StringTokenizer(line, " ");
            st.nextToken();
            String label = st.nextToken();

            if(!hasCord)
            {
                x = -1;//(int)(Math.random() * 500D) + 20;
                y = -1;//(int)(Math.random() * 500D) + 20;
            } else
            {
                String cord = st.nextToken();
                x = Double.parseDouble(cord.substring(2));
                cord = st.nextToken();
                y = Integer.parseInt(cord.substring(2));
            }
            Node node = new Node(label, nIndex);
            if (x!=-1 && y!=-1) {
            	node.setX(x);
            	node.setY(y);
            }
            graph.nodes.add(node);
            nIndex++;
        }

        for(String line = reader.readLine(); line != null; line = reader.readLine())
        {
            if(line.startsWith("*"))
                break;
            StringTokenizer st = new StringTokenizer(line, "-");
            int n1 = Integer.parseInt(st.nextToken()) - 1;
            int n2 = Integer.parseInt(st.nextToken()) - 1;
            Edge edge = new Edge((Node)(graph.nodes.get(n1)), (Node)(graph.nodes.get(n2)), isOrdered);
            graph.edges.add(edge);
        }
        
        System.out.println(filename + " " + graph.getNodesNum() + " --- " + graph.getEdgesNum());
        reader.close();
        return graph;
	}
	
	public void writeGraph(File file) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(file));
        for(int i = 0; i < nodes.size(); i++)
        {
            bw.write((i + 1) + ") " + (Node)nodes.elementAt(i));
            bw.newLine();
        }

        bw.write("edge");
        bw.newLine();
        for(int i = 0; i < edges.size(); i++)
        {
            Edge edge = (Edge)edges.elementAt(i);
            int s = nodes.indexOf(edge.getStartNode()) + 1;
            int e = nodes.indexOf(edge.getEndNode()) + 1;
            bw.write(s + "-" + e);
            bw.newLine();
        }

        bw.write("*");
        bw.newLine();
        bw.close();
	}

    public int getNodesNum(){
    	return this.nodes.size();
    }
    
    public Node getNode(int index) throws Exception{
    	return (Node) (this.nodes.elementAt(index));
    }
    
    public int getEdgesNum(){
    	return this.edges.size();
    }
    
    public Edge getEdge(int index) throws Exception{
    	return (Edge) (this.edges.elementAt(index));
    }
    
    public int getType(){
    	return this.type;
    }
    
    public void setType(int type){
    	this.type = type;
    }    
}
