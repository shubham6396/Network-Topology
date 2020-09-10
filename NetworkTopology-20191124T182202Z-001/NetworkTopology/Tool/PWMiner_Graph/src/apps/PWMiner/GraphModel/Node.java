package apps.PWMiner.GraphModel;

/**
 * @author ITCS 6114
 * Node Model
 */


public class Node{
	int index;
    String label;
    String title;   
    
    int type; // Waiting for using

    double x = -1;
    double y = -1;
    
    public Node(String label, int index){
        this.label = label;
        this.title = label;
        this.index = index;
    }

    public String getLabel(){
        return label;
    }
    
    public String getTitle(){
        return title;
    }

    public int getIndex(){
        return index;
    }
    
    public String toString(){
        return label  + ((x!=-1 && y!=-1)? " x=" + x + " y=" + y : "");
    }

    public double getX(){
        return x;
    }

    public double getY(){
        return y;
    }
    
    public void setX(double x){
    	this.x = x;
    }
    
    public void setY(double y){
    	this.y = y;
    }
    
    public int getType(){
    	return this.type;
    }
    
    public void setType(int type){
    	this.type = type;
    }
}
