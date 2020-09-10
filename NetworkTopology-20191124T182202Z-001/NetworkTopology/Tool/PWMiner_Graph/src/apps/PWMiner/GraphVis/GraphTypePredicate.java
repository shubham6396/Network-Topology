package apps.PWMiner.GraphVis;

import prefuse.data.Schema;
import prefuse.data.Tuple;
import prefuse.visual.VisualItem;
import prefuse.visual.expression.InGroupPredicate;

public class GraphTypePredicate extends InGroupPredicate {
    
    /**
     * Create a new GraphTypePredicate.
     */
    public GraphTypePredicate() {
    }
    
    /**
     * Create a new GraphTypePredicate.
     * @param group @param group the data group name to use as a parameter
     */
    public GraphTypePredicate(int type) {
        super("graph");
    }
    
    /**
     * @see prefuse.data.expression.Expression#get(prefuse.data.Tuple)
     */
    public Object get(Tuple t) {
        return getBoolean(t) ? Boolean.TRUE : Boolean.FALSE;
    }
    
    /**
     * @see prefuse.data.expression.Expression#getBoolean(prefuse.data.Tuple)
     */
    public boolean getBoolean(Tuple t) {
        if ( !(t instanceof VisualItem) )
            return false;
        
        String group = getGroup(t);
        if ( group == null ) {
            return false;
        }
        VisualItem item = (VisualItem)t;
        return item.getVisualization().isInGroup(item, group);
    }

    /**
     * @see prefuse.data.expression.Function#getName()
     */
    public String getName() {
        return "INGROUP";
    }

    /**
     * @see prefuse.data.expression.Expression#getType(prefuse.data.Schema)
     */
    public Class getType(Schema s) {
        return boolean.class;
    }
}
