package apps.PWMiner.GraphVis;

import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.FileOutputStream;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JSplitPane;

import apps.PWMiner.common.Define;

import prefuse.controls.DragControl;
import prefuse.controls.PanControl;
import prefuse.controls.WheelZoomControl;
import prefuse.controls.ZoomControl;
import prefuse.controls.ZoomToFitControl;
import prefuse.util.ui.JPrefuseApplet;
import prefuse.util.ui.UILib;

public class PWGraphHandleApplet extends JPrefuseApplet {
	public void init(){
        UILib.setPlatformLookAndFeel();
        JComponent graphview = buildContentPane();
        this.getContentPane().add(graphview);

        this.setSize(getVisualization().getWidth()+5,getVisualization().getHeight()+50);    	
        
        this.getContentPane().setLayout(new BoxLayout(this.getContentPane(), BoxLayout.Y_AXIS));       
    }
	
    //protected static Logger logger = Logger.getLogger(GraphVisHandler.class);

    private static boolean isLaunched = false;
    static PWGraphHandleApplet imageLib = new PWGraphHandleApplet();
    
    private static final String hover = "hover";   
    public static final String GRAPH = "graph";   
    private static final String nodes = "graph.nodes";
    
	public PWGraphHandleApplet() {
		// TODO Auto-generated constructor stub
	}
	
	public static void isLaunched() {
        if (isLaunched) {
            return;
        }
        isLaunched = true;
    }
	
	public static PWGraphSearchViz getImage(String basedir, String graphFile, boolean isOrdered, String imageName, int imageL, int imageW){
		PWGraphSearchViz viz=null;
		isLaunched();
    	try{
    		viz = new PWGraphSearchViz(basedir + "/" + Define.DATA_DIR + "/" + graphFile, isOrdered, imageL, imageW, 0);
    	    //ZoomToFitAction fitaction = new ZoomToFitAction(viz.getVisualization());
    	    //fitaction.run(1.0);
    		System.out.println("output " + (viz==null? "T" : "F"));
    		if (viz.isVisible()){
	    	    while ( ! viz.isSaveImageEnable() ){
	    	    	
	    	    }
	    	    	
	    	    System.out.println("output = " + basedir  +"/" + Define.GRAPH_VIZ_DIR + "/" + imageName + Define.VIZ_G_SUFFIX);
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
		isLaunched();
    	try{
    		viz = new PWGraphSearchViz(basedir + "/" + Define.DATA_DIR + "/" + graphFile, isOrdered, imageL, imageW, type);
     		if (! viz.isVisible())
    			viz=null;
	    }catch (Exception ex){
	    	ex.printStackTrace();
	    }     
	    return viz;
    }	
	
	public static void saveImageByType(PWGraphSearchViz viz, String basedir, String graphFile, boolean isOrdered, String imageName, int imageL, int imageW, int type){
		isLaunched();
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
	
	static PWGraphSearchViz viz;
	private PWGraphSearchViz getVisualization(){
		return this.viz;
	}
    public static JSplitPane buildContentPane() {
    	//System.setProperty("java.awt.headless","true");
    	
    	//ConfigFile.init();
    	//String pw = "Escherichia_coli_K12/superpathway_of_threonine_metabolism";
/*
    	final String pw = "Escherichia_coli_K12/b";
		final String graphFile = pw + ".grp";
    	viz = GraphVisHandler.getImageByType(System.getProperty("user.dir") + "/..", graphFile, true, pw, 400, 400, CommonDef.ISMAPPING);
*/    	
/*
    	final String pw = "Escherichia_coli_K12/b1";
		final String graphFile = pw + ".grp";
    	viz = GraphVisHandler.getImageByType(System.getProperty("user.dir") + "/..", graphFile, true, pw, 400, 400, CommonDef.ISPATTERN);
*/    	 	
    	final String pw = "Escherichia_coli_K12/superpathway_of_leucine-_valine-_and_isoleucine_biosynthesis";//b3";
		final String graphFile = pw + ".grp";
    	viz = GraphVisHandler.getImageByType(System.getProperty("user.dir") + "/..", graphFile, true, pw, 400, 400, CommonDef.ISTEXT); 
    	
    	//JFrame frame = new JFrame("Graph Visualization");
        
        //main display controls
        viz.setSize(400,400);
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
				GraphVisHandler.saveImageByType(viz, System.getProperty("user.dir") + "/..", graphFile, true, pw, 400, 400, CommonDef.ISMAPPING);
				System.out.println("---------");
				
				//Upload
				
			}
    		
    	});
    	fpanel.add(saveButton);
    	
    	JPanel fpane2 = new JPanel();
    	fpane2.add(new JButton("Try"));
    	
    	JSplitPane split = new JSplitPane(JSplitPane.VERTICAL_SPLIT, viz, fpanel);
    	split.setOneTouchExpandable(true);
    	split.setContinuousLayout(false);
    	split.setDividerLocation(400);

    	// now we run our action list
    	viz.getVisualization().run("draw");

    	//frame.getContentPane().add(split);
        //frame.setSize(viz.getWidth()+10,viz.getHeight()+160);    	
        
        //frame.pack();
        //frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        
        //<-CQ ADD
        //frame.getContentPane().setLayout(new BoxLayout(frame.getContentPane(), BoxLayout.Y_AXIS));
        // -->
        
        //frame.setVisible(true);
		return split;
        
    }	
}
