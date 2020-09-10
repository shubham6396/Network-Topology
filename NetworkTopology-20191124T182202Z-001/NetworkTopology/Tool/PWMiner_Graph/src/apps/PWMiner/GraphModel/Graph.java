package apps.PWMiner.GraphModel;

import java.io.File;
import java.io.IOException;

/**
 * @author ITCS 6114
 * Graph Model Interface
 */
public interface Graph {
    
    public abstract void writeGraph(File file) throws IOException;
    
    public abstract void loadGraph(Node nodesArr[], Edge edgesArr[]);
    
    public abstract Graph readGraph(String filename, boolean isOrdered) throws IOException;
    
    public abstract int getNodesNum();
    
    public abstract Node getNode(int index) throws Exception;
    
    public abstract int getEdgesNum();
    
    public abstract Edge getEdge(int index) throws Exception;
}
