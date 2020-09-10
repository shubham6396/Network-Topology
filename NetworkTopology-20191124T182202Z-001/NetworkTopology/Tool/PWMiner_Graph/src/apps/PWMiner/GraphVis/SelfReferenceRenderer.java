package apps.PWMiner.GraphVis;

import java.awt.Shape;
import java.awt.geom.Ellipse2D;

import prefuse.render.EdgeRenderer;
import prefuse.visual.EdgeItem;
import prefuse.visual.VisualItem;

public class SelfReferenceRenderer extends EdgeRenderer { 
	private Ellipse2D m_ellipse = new Ellipse2D.Float(); 
	protected Shape getRawShape(VisualItem item) { 
		try { 
			EdgeItem edge = (EdgeItem) item; 
			VisualItem item1 = edge.getSourceItem(); 
			VisualItem item2 = edge.getTargetItem();  
			 
			//self interaction 
			if (item1 == item2) {  
				getAlignedPoint(m_tmpPoints[0], item1.getBounds(), m_xAlign1, m_yAlign1); 
				getAlignedPoint(m_tmpPoints[1], item2.getBounds(), m_xAlign2, m_yAlign2); 
				 
				m_curWidth = (int) Math.round(m_width * getLineWidth(item)); 
				m_ellipse.setFrame(m_tmpPoints[0].getX(), m_tmpPoints[0].getY(), 40, 30);  
				return m_ellipse; 
			} 
		} 
		catch(Exception ex) { 
			ex.printStackTrace(); 
			return null; 
		} 
	 
		return super.getRawShape(item); 
	} //getRawShape 
} 
