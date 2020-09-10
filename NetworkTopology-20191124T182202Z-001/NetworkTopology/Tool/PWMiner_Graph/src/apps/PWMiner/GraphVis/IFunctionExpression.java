package apps.PWMiner.GraphVis;

import prefuse.data.Node;
import prefuse.data.Schema;
import prefuse.data.Tuple;
import prefuse.data.expression.FunctionExpression;

/**
 * FunctionTable.addFunction("IFunctionExpression", IFunctionExpression.class);
 * @author Qiong
 *
 */
public class IFunctionExpression extends FunctionExpression {
	public IFunctionExpression() {
		super(0);
		
	}
	public String getName() {
		// TODO Auto-generated method stub
		return "IFunctionExpression";
	}

	public Class getType(Schema s) {
        return int.class;//Example
    }
	
	public int getInt(Tuple t) { //Example
        return (t instanceof Node ? ((Node)t).getInDegree() : 0 );
    }
}
