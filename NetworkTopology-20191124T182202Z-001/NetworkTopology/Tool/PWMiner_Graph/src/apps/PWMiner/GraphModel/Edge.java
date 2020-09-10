package apps.PWMiner.GraphModel;

/**
 * @author ITCS6114
 * Edge Model
 */
public class Edge{
    Node startNode;
    Node endNode;
    boolean isOrdered;

    public Edge(Node startNode, Node endNode, boolean isOrdered){
        this.startNode = startNode;
        this.endNode = endNode;
        this.isOrdered = isOrdered;
    }

    public Node getStartNode(){
        return startNode;
    }

    public Node getEndNode(){
        return endNode;
    }

    public String toString(){
        return " < " + startNode + "," + endNode + "," + isOrdered + " > ";
    }
}
