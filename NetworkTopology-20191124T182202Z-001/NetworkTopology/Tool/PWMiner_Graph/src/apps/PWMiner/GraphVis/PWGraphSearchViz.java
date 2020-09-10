package apps.PWMiner.GraphVis;

import java.awt.Color;
import java.awt.Cursor;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.geom.Rectangle2D;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterException;
import java.awt.print.PrinterJob;
import java.io.InputStream;
import java.util.Iterator;

import javax.swing.AbstractAction;
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
import prefuse.data.Tuple;

import prefuse.data.event.TupleSetListener;
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
import prefuse.visual.expression.InGroupPredicate;

public class PWGraphSearchViz extends Display implements ActivityListener,Printable {
    //protected static Logger logger = Logger.getLogger(PWGraphSearchViz.class);	
    private int arrowHeadWidth = 8;
    private int arrowHeadHeight = 6;
    private boolean directedGraph = true;

    //private int edgeType = prefuse.Constants.EDGE_TYPE_CURVE;
	private int edgeType = prefuse.Constants.EDGE_TYPE_LINE;

    private int arrowType = prefuse.Constants.EDGE_ARROW_FORWARD;
    //private int arrowType = prefuse.Constants.EDGE_ARROW_REVERSE;
    
    public static final String GRAPH = "graph";
    private static final String hover = "hover";    
    public static final String NODES = "graph.nodes";
    public static final String EDGES = "graph.edges";
    private static final String linear = "linear";
    private static final String SIZE = "size";
    
    private boolean enableSaveImage = true;
    
    private int type = CommonDef.ISNULL;
    //private Visualization m_vis = null;
    // MAD - the name of the different neighborgroups
    private static final String[] neighborGroups = { "sourceNode", "targetNode", "bothNode" };
    
    public PWGraphSearchViz(String GraphFile, boolean isOrdered, int length, int width, int type) throws Exception{
        super(new Visualization());
        
        // generate a graph
        initGraph(GraphFile, isOrdered);
        
        buildViz(length, width, type);
    }
    
    public PWGraphSearchViz(InputStream in, boolean isOrdered, int length, int width, int type) throws Exception{
        super(new Visualization());
        
        // generate a graph
        initGraphByInputStream(in, isOrdered);
        
        buildViz(length, width, type);
    }
    
    private void buildViz(int length, int width, int type) throws Exception{
        this.setSize(length, width);
        //this.setLocation(5, 10);
        this.setBounds(5,5, width-10, width-10);        
        this.type = type;
        //System.out.println("--- width=" + this.getWidth() + " height=" + this.getHeight() + " graphfile=" + GraphFile);
        this.m_vis.setValue(EDGES, null, VisualItem.INTERACTIVE, Boolean.TRUE);
        
    	// <- CQ Add 01-03-2008 fix selected focus nodes
    	TupleSet focusGroup = m_vis.getGroup(Visualization.FOCUS_ITEMS); 
    	focusGroup.addTupleSetListener(new TupleSetListener() {
    	    public void tupleSetChanged(TupleSet ts, Tuple[] add, Tuple[] rem){
    			for ( int i=0; i<rem.length; ++i )
    			    ((VisualItem)rem[i]).setFixed(false);
    			for ( int i=0; i<add.length; ++i ) {
    			    ((VisualItem)add[i]).setFixed(false);
    			    ((VisualItem)add[i]).setFixed(true);
    			}
    			if ( ts.getTupleCount() == 0 ) {
    			    ts.addTuple(rem[0]);
    			    ((VisualItem)rem[0]).setFixed(false);
    			}
    			//m_vis.run("draw");//CQ 01-23
    	    }
    	});

    	// --------------------------------------------------------------------
    	// create actions to process the visual data  	
    	// CQ ->
    	
        //m_vis.getVisualItem(NODES, m_vis.g).setFixed(true);
        TupleSet ts=m_vis.getVisualGroup(NODES);
        for ( Iterator items = ts.tuples(); items.hasNext(); ) {
            VisualItem item = (VisualItem)items.next();
            item.setInt(CommonDef.TYPE, this.type);
            if (this.type==CommonDef.ISMAPPING){
            	if (item.getString(CommonDef.LABEL).indexOf("\n") == -1)
            		item.setInt(CommonDef.TYPE, CommonDef.ISTEXT);
            }
            item.setInt(CommonDef.PATTERN, ColorLib.rgb(200,200,255));
            item.setInt(CommonDef.TEXT, ColorLib.gray(255));
            //item.setFixed(true);
        }
        // standard labelRenderer for the given label
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
        
        // <- CQ Add 01-03-2008
        nodeFill.add(VisualItem.FIXED, ColorLib.rgb(255,100,100));
        nodeFill.add(VisualItem.HIGHLIGHT, ColorLib.rgb(255,200,125));
        //nodeFill.add(CommonDef.PATTERN, ColorLib.rgb(200,200,255));
        //nodeFill.add(CommonDef.TEXT, ColorLib.gray(255));
        
        nodeFill.add("ingroup('_search_')", ColorLib.rgb(255,190,190));        
        //MAD - here we define the colors in which the neighbour nodes should be filled
        nodeFill.add(new InGroupPredicate(neighborGroups[0]), ColorLib.rgb(0, 0,
    			250));
        nodeFill.add(new InGroupPredicate(neighborGroups[1]), ColorLib.rgb(0, 250,
    			0));
        nodeFill.add(new InGroupPredicate(neighborGroups[2]), ColorLib.rgb(250, 0,
    			0));
    	// CQ ->
        
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
    	// <- CQ ADD 01-03-2008
        //m_vis.putAction("onlydraw", draw); //draw only
        
    	// CQ ->
        
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
        flayout.setLayoutBounds(new Rectangle2D.Double(5, 5, width-5, width-5));
        // <- CQ ADD 01-03-2008
        flayout.setIterations(100);
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
        //this.m_vis.run("draw");
//ForceSimulator fsim = flayout.getForceSimulator();
//fsim.addForce(new SpringForce(-(2E-2f), 100f));
	}
	
    private void initGraph(String GraphFile, boolean isOrdered) throws Exception {
    	//ReadPWPlainTextGraph
    	PWPlainGraph pwGraph = PWPlainGraph.readGraph(GraphFile, isOrdered);
    	buildGraph(pwGraph, isOrdered);
    }
    
    private void initGraphByInputStream(InputStream in, boolean isOrdered) throws Exception {
    	//ReadPWPlainTextGraph
    	PWPlainGraph pwGraph = PWPlainGraph.readGraphByInputStream(in, isOrdered);
    	buildGraph(pwGraph, isOrdered);
    }
    
    private void buildGraph(PWPlainGraph pwGraph, boolean isOrdered) throws Exception {
    	pwGraph.setType(this.type);
    	
    	//build prefuse.data.graph
    	prefuse.data.Graph graph = new prefuse.data.Graph();
        graph.setDirected(isOrdered);
        // first define columns
        graph.addColumn(CommonDef.LABEL, String.class);
        graph.addColumn(CommonDef.TYPE, int.class);
        graph.addColumn(CommonDef.PATTERN, int.class);
        graph.addColumn(CommonDef.TEXT, int.class);
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
            
            //System.out.println(nod.getLabel() + " \n");
            if ((x>0) && (x<nlab.length()-1)) h.setString(CommonDef.LABEL, nod.getLabel().replaceAll("=", "\n"));
            else if ((x==0) || (x==nlab.length()-1))
            	h.setString(CommonDef.LABEL, nod.getLabel().replaceAll("=", ""));
            else
            	h.setString(CommonDef.LABEL, nod.getLabel());
            
            if (pwGraph.getType() == CommonDef.ISPATTERN)
            	h.setInt(CommonDef.TYPE, CommonDef.ISPATTERN);
            else if (pwGraph.getType() == CommonDef.ISTEXT)
            	h.setInt(CommonDef.TYPE, CommonDef.ISTEXT);
            else if (pwGraph.getType() == CommonDef.ISMAPPING){
            	if (nod.getLabel().indexOf("\n") == -1)
            		h.setInt(CommonDef.TYPE, CommonDef.ISTEXT);
            	else
            		h.setInt(CommonDef.TYPE, CommonDef.ISMAPPING);
            	
            }
            //System.out.println(h.getString(LABEL) + " \n");
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

    /**
     * Print the current graph.
     */
    class PrintAction extends AbstractAction {
    	public static final String printKey = "print";
    	PrintAction() {
    		super(printKey);
    	}

    	public void actionPerformed(ActionEvent e) {
    		PrinterJob printJob = PrinterJob.getPrinterJob();
    		printJob.setPrintable(PWGraphSearchViz.this );
    		if (printJob.printDialog()) {
    			try {
    				//boolean oldvalue = false;
    				//double oldscale = graph.getScale();
    				//setScale(1);
    				//oldvalue = graph.isPageVisible();
    				printJob.print();
    				//setScale(oldscale);
    				//graph.setPageVisible(oldvalue);
    			} catch (Exception printException) {
    				printException.printStackTrace();
    			}
    		}
    	}
    }

	public int print(Graphics arg0, PageFormat arg1, int arg2) throws PrinterException {
		// TODO Auto-generated method stub
		return 0;
	}
}
