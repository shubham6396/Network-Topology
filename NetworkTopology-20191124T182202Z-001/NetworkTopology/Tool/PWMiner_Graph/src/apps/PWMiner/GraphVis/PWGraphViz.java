package apps.PWMiner.GraphVis;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Font;
import java.awt.geom.Rectangle2D;
import java.util.Iterator;
import apps.PWMiner.GraphModel.PWPlainGraph;
import prefuse.Constants;
import prefuse.Display;
import prefuse.Visualization;
import prefuse.action.ActionList;
import prefuse.action.RepaintAction;
import prefuse.action.assignment.ColorAction;
import prefuse.action.assignment.FontAction;
import prefuse.activity.Activity;
import prefuse.activity.ActivityListener;
import prefuse.data.Node;
import prefuse.data.expression.Predicate;
import prefuse.data.expression.parser.ExpressionParser;
import prefuse.data.tuple.TupleSet;
import prefuse.render.AbstractShapeRenderer;
import prefuse.render.DefaultRendererFactory;
import prefuse.render.EdgeRenderer;
import prefuse.render.LabelRenderer;
import prefuse.util.ColorLib;
import prefuse.util.GraphicsLib;
import prefuse.util.display.DisplayLib;
import prefuse.util.display.ItemBoundsListener;
import prefuse.visual.VisualItem;

public class PWGraphViz extends Display implements ActivityListener {
    //protected static Logger logger = Logger.getLogger(PWGraphViz.class);	
    private int arrowHeadWidth = 8;
    private int arrowHeadHeight = 6;
    private boolean directedGraph = true;

    //private int edgeType = prefuse.Constants.EDGE_TYPE_CURVE;
	private int edgeType = prefuse.Constants.EDGE_TYPE_LINE;

    private int arrowType = prefuse.Constants.EDGE_ARROW_FORWARD;
    //private int arrowType = prefuse.Constants.EDGE_ARROW_REVERSE;
    
    private static final String LABEL = "label";

    public static final String GRAPH = "graph";
    private static final String hover = "hover";    
    public static final String NODES = "graph.nodes";
    public static final String EDGES = "graph.edges";
    private static final String linear = "linear";
    private static final String SIZE = "size";
    
    private boolean enableSaveImage = true;
    
    // MAD - the name of the different neighborgroups
    private static final String[] neighborGroups = { "sourceNode", "targetNode", "bothNode" };
    
    public PWGraphViz(String GraphFile, boolean isOrdered, int length, int width) throws Exception{
        super(new Visualization());
        
        this.setSize(length, width);
        System.out.println("--- width=" + this.getWidth() + " height=" + this.getHeight());
        this.m_vis.setValue(EDGES, null, VisualItem.INTERACTIVE, Boolean.TRUE);
        
        // generate a graph
        initGraph(GraphFile, isOrdered);
    	
        //m_vis.getVisualItem(NODES, m_vis.g).setFixed(true);
        TupleSet ts=m_vis.getVisualGroup(NODES);
        for ( Iterator items = ts.tuples(); items.hasNext(); ) {
            VisualItem item = (VisualItem)items.next();
            //item.setFixed(true);
        }
        // standard labelRenderer for the given label
        //LabelRenderer nodeRenderer = new LabelRenderer(LABEL);//CQ 01-23
        ILabelRenderer nodeRenderer = new ILabelRenderer(CommonDef.LABEL);        
        nodeRenderer.setRenderType(AbstractShapeRenderer.RENDER_TYPE_DRAW_AND_FILL); 
        nodeRenderer.setHorizontalAlignment(Constants.LEFT); 
        nodeRenderer.setRoundedCorner(8,8);         
    	EdgeRenderer edgeRenderer = new SelfReferenceRenderer();
    	edgeRenderer.setArrowHeadSize(arrowHeadWidth, arrowHeadHeight);

    	edgeRenderer.setEdgeType(edgeType); 
    	edgeRenderer.setArrowType(arrowType);

        // rendererFactory for the visualization items 
        DefaultRendererFactory rendererFactory = new DefaultRendererFactory();
        
        // set the labelRenderer
        rendererFactory.setDefaultRenderer(nodeRenderer);
        rendererFactory.setDefaultEdgeRenderer(edgeRenderer);        
        this.m_vis.setRendererFactory(rendererFactory);

        // Color Actions
        ColorAction nodeText = new ColorAction(NODES, VisualItem.TEXTCOLOR);
        nodeText.setDefaultColor(ColorLib.gray(0));
        ColorAction nodeStroke = new ColorAction(NODES, VisualItem.STROKECOLOR);
        nodeStroke.setDefaultColor(ColorLib.gray(100));
        ColorAction nodeFill = new ColorAction(NODES, VisualItem.FILLCOLOR);
        nodeFill.setDefaultColor(ColorLib.gray(255));
        ColorAction edgeStrokes = new ColorAction(EDGES, VisualItem.STROKECOLOR);
        edgeStrokes.setDefaultColor(ColorLib.gray(100));      
            
        //DataSizeAction
        //DataSizeAction nodeDataSizeAction = new DataSizeAction(NODES, SIZE);
        //nodeDataSizeAction.setIs2DArea(true);
        
        // bundle the color actions
        ActionList draw = new ActionList();
        draw.add(nodeText);
        draw.add(nodeStroke);
        draw.add(nodeFill);
        draw.add(edgeStrokes);
        draw.add(new ColorAction(EDGES, VisualItem.FILLCOLOR, ColorLib.gray(100)));   
        //draw.add(nodeDataSizeAction);
        
        //draw graph root with a different color
        Predicate rootPredict =  (Predicate)ExpressionParser.parse("INDEGREE()=0");//treedepth()=1"); 
        ColorAction rootNodeColor =new ColorAction("graph.nodes",rootPredict,VisualItem.FILLCOLOR, ColorLib.color(Color.YELLOW)); 
        draw.add(rootNodeColor);
        
        //draw nodes font
        FontAction nodeFont =new FontAction("graph.nodes", new Font("Serif", Font.PLAIN, 14) );
        draw.add(nodeFont);
        
        this.m_vis.putAction("draw", draw);
        
        // create the layout action for the graph // run actions
        //NodeLinkTreeLayout treeLayout = new NodeLinkTreeLayout(GRAPH); //CQ
        //this.m_vis.putAction("treeLayout", treeLayout); //CQ
        
        //this.m_vis.run("treeLayout"); //CQ
        
        /*DotLayout dotLayout = new DotLayout();
        this.m_vis.putAction("Layout", dotLayout);
        this.m_vis.run("layout");*/
        
        ActionList animate = new ActionList(8000);//Activity.INFINITY);
        //PolarLocationAnimator animate = new PolarLocationAnimator(GRAPH);
        
        IForceDirectedLayout flayout= new IForceDirectedLayout(GRAPH, true, false);
        // <- CQ ADD 01-03-2008
        //flayout.setIterations(100);
        // CQ ->
        //flayout.setLayoutBounds(new Rectangle(0, 0, 400, 400));
        animate.add(flayout);
        animate.add(new RepaintAction());
        animate.addActivityListener(this);
        m_vis.putAction("layout", animate);

        /*
        RadialTreeLayout treeLayout = new RadialTreeLayout(tree);
        //treeLayout.setAngularBounds(-Math.PI/2, Math.PI);
        m_vis.putAction("treeLayout", treeLayout);
        
        CollapsedSubtreeLayout subLayout = new CollapsedSubtreeLayout(tree);
        m_vis.putAction("subLayout", subLayout);
        
        // create the filtering and layout
        ActionList filter = new ActionList();
        filter.add(new TreeRootAction(tree));
        filter.add(fonts);
        filter.add(treeLayout);
        filter.add(subLayout);
        filter.add(textColor);
        filter.add(nodeColor);
        filter.add(edgeColor);
        m_vis.putAction("filter", filter);
        
        // animated transition
        ActionList animate = new ActionList(1250);
        animate.setPacingFunction(new SlowInSlowOutPacer());
        animate.add(new QualityControlAnimator());
        animate.add(new VisibilityAnimator(tree));
        animate.add(new PolarLocationAnimator(treeNodes, linear));
        animate.add(new ColorAnimator(treeNodes));
        animate.add(new RepaintAction());
        m_vis.putAction("animate", animate);
        m_vis.alwaysRunAfter("filter", "animate");
         */
        
        m_vis.run("layout");  //runAfter("draw", "layout");
        this.m_vis.run("draw");
        //ForceSimulator fsim = flayout.getForceSimulator();
        //fsim.addForce(new SpringForce(-(2E-2f), 100f));
        
        
	}
	
    private void initGraph(String GraphFile, boolean isOrdered) throws Exception {
    	//ReadPWPlainTextGraph
    	PWPlainGraph pwGraph = PWPlainGraph.readGraph(GraphFile, isOrdered);
    	
    	//build prefuse.data.graph
    	prefuse.data.Graph graph = new prefuse.data.Graph();
        graph.setDirected(isOrdered);
        // first define a column
        graph.addColumn(LABEL, String.class);
        // then add nodes
        
        if (pwGraph.getNodesNum()==0)
        	isVisible=false;
        for(int i = 0; i < pwGraph.getNodesNum(); i++){
            apps.PWMiner.GraphModel.Node nod = pwGraph.getNode(i);
            Node h = graph.addNode();
            String nlab = nod.getLabel();
            int x = nlab.indexOf('=');
            /*
            int len = (x > nlab.length() -x-1 ? x : nlab.length() -x-1);
            char[] padding = new char[len+2];
            for (int j=0; j<len+2; j++)
            	padding[j]='-';
            h.setString(LABEL, nod.getLabel().replaceAll("=", "\n"+ new String(padding) + "\n"));
            */
            
            System.out.println(nod.getLabel() + " \n");
            if (x>-1) h.setString(LABEL, nod.getLabel().replaceAll("=", "\n"));
            else h.setString(LABEL, nod.getLabel());
            System.out.println(h.getString(LABEL) + " \n");
        }        
        for(int j = 0; j < pwGraph.getEdgesNum(); j++){
            apps.PWMiner.GraphModel.Edge edg = pwGraph.getEdge(j);
            graph.addEdge(graph.getNode(edg.getStartNode().getIndex()), graph.getNode(edg.getEndNode().getIndex()));
        }    

        this.m_vis.addGraph(GRAPH, graph);

    }
    
    private boolean isVisible=true;
    public boolean isVisible(){
    	return this.isVisible;
    }
    

	
	/**
	 * Updates the view making sure all the colors etc. are correct.
	 */
	public void runFilterUpdate() {
		m_vis.run("update");
	}
	
	/**
	 * Schedules the dot layout to be performed immediately
	 */
	private void runDotLayout() {
	    setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
	    m_vis.run("dotLayout");	    
	} //
	
	/**
	 * Runs the layout. All nodes and edges are moved to the positions determined
	 * by the used layout algorithm. Currently the layout is performed by dot.
	 */
	public void runLayout() {
		runDotLayout();
	}
	
	public boolean isSaveImageEnable(){
		return this.enableSaveImage;
	}
	
	//----Start of implementing ActivityManager----------------------------
    /**
     * Called when an activity has been scheduled with an ActivityManager
     * @param a the scheduled Activity
     */
    public void activityScheduled(Activity a){
    	enableSaveImage = false;
    }
    
    /**
     * Called when an activity is first started.
     * @param a the started Activity
     */
    public void activityStarted(Activity a){

    }
    
    /**
     * Called when an activity is stepped.
     * @param a the stepped Activity
     */
    public void activityStepped(Activity a){
    	
    }
    
    /**
     * Called when an activity finishes.
     * @param a the finished Activity
     */
    public void activityFinished(Activity a){
    	enableSaveImage = true;
    }
    
    /**
     * Called when an activity is cancelled.
     * @param a the cancelled Activity
     */
    public void activityCancelled(Activity a){
    	
    }
	//----End of implementing ActivityManager----------------------------
    
    public class MyLabelRenderer extends LabelRenderer {
  	    public String getText(VisualItem vi) {
  	        String name = vi.getString("name");

  	        if (name.indexOf("=") >-1 ){
	  	    	name.replaceAll("=", "\n");
	    	    System.out.println(name);
	    	    vi.setString("name", name.replaceAll("=", "\n"));
  	        }
  	        return name + "\n" + vi.getString("gender");
  	    }
    }
    
    public static class FitOverviewListener implements ItemBoundsListener {
    	private Rectangle2D m_bounds = new Rectangle2D.Double();
    	private Rectangle2D m_temp = new Rectangle2D.Double();
    	private double m_d = 15;
    	public void itemBoundsChanged(Display d) {
    	    d.getItemBounds(m_temp);
    	    GraphicsLib.expand(m_temp, 25/d.getScale());

    	    double dd = m_d/d.getScale();
    	    double xd = Math.abs(m_temp.getMinX()-m_bounds.getMinX());
    	    double yd = Math.abs(m_temp.getMinY()-m_bounds.getMinY());
    	    double wd = Math.abs(m_temp.getWidth()-m_bounds.getWidth());
    	    double hd = Math.abs(m_temp.getHeight()-m_bounds.getHeight());
    	    if ( xd>dd || yd>dd || wd>dd || hd>dd ) {
    		m_bounds.setFrame(m_temp);
    		DisplayLib.fitViewToBounds(d, m_bounds, 0);
    	    }
    	}
    }

}
