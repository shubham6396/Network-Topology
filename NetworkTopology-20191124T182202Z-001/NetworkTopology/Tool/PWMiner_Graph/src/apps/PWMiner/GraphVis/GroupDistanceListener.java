package apps.PWMiner.GraphVis;

import java.util.ArrayList;
import java.util.Iterator;

import prefuse.Display;
import prefuse.Visualization;
import prefuse.data.Tuple;
import prefuse.data.event.TupleSetListener;
import prefuse.data.tuple.TupleSet;
import prefuse.visual.VisualItem;

/**
 *	Generic TupleSetListener applying a distance filter to a given focus group 
 *  If the focus group is empty, the filter is disabled. (you might modify this 
 *  according to your needs)
 *
 * @author martin dudek
 *
 */

public class GroupDistanceListener implements TupleSetListener {

	String graph;
	Visualization vis;
	GraphDistanceFilter filter;
	String drawAction;
	
	ArrayList previousVisibleItems;
	int distance;
	
	boolean lastTimeFiltered = false;
	Display display = null;
	
	public GroupDistanceListener(String graphName,Display display,GraphDistanceFilter filter,String drawAction) {
	    this(graphName,display.getVisualization(),filter,drawAction,1);
	    this.display = display;
	}
	
	public GroupDistanceListener(String graph,Visualization vis,GraphDistanceFilter filter,String drawAction,int distance) {
	    this.graph = graph;
	    this.vis = vis;
	    this.filter = filter;
	    this.drawAction = drawAction;
	    this.distance = distance;
	    previousVisibleItems = new ArrayList();
	}
	
	public void setDistance(int distance) {
	    this.distance = distance;
	    filter.setDistance(distance);
	}
	
	public void tupleSetChanged(TupleSet ts, Tuple[] add, Tuple[] rem){
	    if ( ts.getTupleCount() == 0 ) {
		if (previousVisibleItems != null) {
		    Iterator iter = previousVisibleItems.iterator(); // reconstructimg the pre filtered state
		    while (iter.hasNext()) {
			VisualItem aItem = (VisualItem) iter.next();
			aItem.setVisible(true);
		    }
		}
		lastTimeFiltered = false;
		filter.setEnabled(false);
	
	    } else {
		if (!lastTimeFiltered) { // remembering the last unfiltered set of visible items
		    previousVisibleItems.clear();
		    Iterator iter  = vis.visibleItems(graph);
		    while (iter.hasNext()) {
				VisualItem aItem = (VisualItem) iter.next();
				previousVisibleItems.add(aItem);
		    }
		}
		lastTimeFiltered = true;
		filter.setEnabled(true);
		filter.setDistance(distance);
	
	    }
	    vis.run(drawAction);
	    display.repaint();
	}
}
