package apps.PWMiner.GraphVis;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.geom.RectangularShape;
import java.awt.geom.RoundRectangle2D;

import prefuse.render.AbstractShapeRenderer;
import prefuse.util.ColorLib;
import prefuse.util.GraphicsLib;
import prefuse.visual.VisualItem;

public class IGraphicsLib extends GraphicsLib{

    /**
     * Render a shape associated with a VisualItem into a graphics context. This
     * method uses the {@link java.awt.Graphics} interface methods when it can,
     * as opposed to the {@link java.awt.Graphics2D} methods such as
     * {@link java.awt.Graphics2D#draw(java.awt.Shape)} and
     * {@link java.awt.Graphics2D#fill(java.awt.Shape)}, resulting in a
     * significant performance increase on the Windows platform, particularly
     * for rectangle and line drawing calls.
     * @param g the graphics context to render to
     * @param item the item being represented by the shape, this instance is
     * used to get the correct color values for the drawing
     * @param shape the shape to render
     * @param stroke the stroke type to use for drawing the object.
     * @param type the rendering type indicating if the shape should be drawn,
     * filled, or both. One of
     * {@link prefuse.render.AbstractShapeRenderer#RENDER_TYPE_DRAW},
     * {@link prefuse.render.AbstractShapeRenderer#RENDER_TYPE_FILL},
     * {@link prefuse.render.AbstractShapeRenderer#RENDER_TYPE_DRAW_AND_FILL}, or
     * {@link prefuse.render.AbstractShapeRenderer#RENDER_TYPE_NONE}.
     */
    public static void paint(Graphics2D g, VisualItem item,
                             Shape shape, BasicStroke stroke, int type, int nodeType)
    {
        // if render type is NONE, then there is nothing to do
        if ( type == AbstractShapeRenderer.RENDER_TYPE_NONE )
            return;
        
        // set up colors
        
        Color strokeColor = ColorLib.getColor(item.getStrokeColor());
        Color fillColor = ColorLib.getColor(item.getFillColor());
        Color patternC = ColorLib.getColor(item.getInt(CommonDef.PATTERN));
        Color textC = ColorLib.getColor(item.getInt(CommonDef.TEXT));
   
        boolean sdraw = (type == AbstractShapeRenderer.RENDER_TYPE_DRAW ||
                         type == AbstractShapeRenderer.RENDER_TYPE_DRAW_AND_FILL) &&
                        strokeColor.getAlpha() != 0;
        boolean fdraw = (type == AbstractShapeRenderer.RENDER_TYPE_FILL ||
                         type == AbstractShapeRenderer.RENDER_TYPE_DRAW_AND_FILL) &&
                        fillColor.getAlpha() != 0;
        if ( !(sdraw || fdraw) ) return;
        
        Stroke origStroke = null;
        if ( sdraw ) {
            origStroke = g.getStroke();
            g.setStroke(stroke);
        }
        
        int x, y, w, h, aw, ah;
        double xx, yy, ww, hh;

        // see if an optimized (non-shape) rendering call is available for us
        // these can speed things up significantly on the windows JRE
        // it is stupid we have to do this, but we do what we must
        // if we are zoomed in, we have no choice but to use
        // full precision rendering methods.
        AffineTransform at = g.getTransform();
        double scale = Math.max(at.getScaleX(), at.getScaleY());
        if ( scale > 1.5 ) {
            if (fdraw) { g.setPaint(fillColor);   g.fill(shape); }
            if (sdraw) { g.setPaint(strokeColor); g.draw(shape); }
        }
        else if ( shape instanceof RectangularShape )
        {
            RectangularShape r = (RectangularShape)shape;
            xx = r.getX(); ww = r.getWidth(); 
            yy = r.getY(); hh = r.getHeight();
            
            x = (int)xx;
            y = (int)yy;
            w = (int)(ww+xx-x);
            h = (int)(hh+yy-y);        	
        	if ((item.getInt(CommonDef.TYPE) == CommonDef.ISPATTERN) || (item.getInt(CommonDef.TYPE) == CommonDef.ISTEXT)){
        		//System.out.println(item.getString("label") + " : " + item.getInt("type"));
        		fillColor = ( (item.getInt(CommonDef.TYPE) == CommonDef.ISPATTERN)? patternC : textC);
        		if ( shape instanceof Rectangle2D ) {
                    if (fdraw) {
                        g.setPaint(fillColor);
                        g.fillRect(x, y, w, h);
                    }
                    if (sdraw) {
                        g.setPaint(strokeColor);
                        g.drawRect(x, y, w, h);
                    }
                } else if ( shape instanceof RoundRectangle2D ) {
                    RoundRectangle2D rr = (RoundRectangle2D)shape;
                    aw = (int)rr.getArcWidth();
                    ah = (int)rr.getArcHeight();
                    if (fdraw) {
                        g.setPaint(fillColor);
                        g.fillRoundRect(x, y, w, h, aw, ah);
                    }
                    if (sdraw) {
                        g.setPaint(strokeColor);
                        g.drawRoundRect(x, y, w, h, aw, ah);
                    }
                } else if ( shape instanceof Ellipse2D ) {
                    if (fdraw) {
                        g.setPaint(fillColor);
                        g.fillOval(x, y, w, h);
                    }
                    if (sdraw) {
                        g.setPaint(strokeColor);
                        g.drawOval(x, y, w, h);
                    }
                } else {
                    if (fdraw) { g.setPaint(fillColor);   g.fill(shape); }
                    if (sdraw) { g.setPaint(strokeColor); g.draw(shape); }
                }
        	}else if (item.getInt(CommonDef.TYPE) == CommonDef.ISMAPPING){
        		xx = r.getX(); ww = r.getWidth(); 
                yy = r.getY(); hh = r.getHeight();
                
                int x1 = (int)xx;
                int y1 = (int)yy;
                w = (int)(ww+xx-x1);
                int h1 = (int)(hh/2+yy-y1);
                
                int x2 = (int)xx;
                int y2 = (int)(yy + h1);
                int h2 = (int)(hh + yy -y2);
                
        		if ( shape instanceof Rectangle2D ) {
                    if (fdraw) {
                        g.setPaint(textC);
                        g.fillRect(x1, y1, w, h1);
                        g.setPaint(patternC);
                        g.fillRect(x2, y2, w, h2);
                    }
                    if (sdraw) {
                        g.setPaint(strokeColor);
                        g.drawRect(x1, y1, w, h1);
                        g.drawRect(x2, y2, w, h2);
                    }
                } else if ( shape instanceof RoundRectangle2D ) {
                    RoundRectangle2D rr = (RoundRectangle2D)shape;
                    aw = (int)rr.getArcWidth();
                    ah = (int)rr.getArcHeight();
                    if (fdraw) {
                        g.setPaint(textC);
                        g.fillRoundRect(x1, y1, w, h1, aw, ah);
                        g.setPaint(patternC);
                        g.fillRoundRect(x2, y2, w, h2, aw, ah);
                    }
                    if (sdraw) {
                        g.setPaint(strokeColor);
                        g.drawRoundRect(x1, y1, w, h1, aw, ah);
                        g.drawRoundRect(x2, y2, w, h2, aw, ah);                        
                    }
                } else if ( shape instanceof Ellipse2D ) {
                    if (fdraw) {
                        g.setPaint(textC);
                        g.fillOval(x1, y1, w, h1);
                        g.setPaint(patternC);
                        g.fillOval(x2, y2, w, h2);                        
                    }
                    if (sdraw) {
                        g.setPaint(strokeColor);
                        g.drawOval(x1, y1, w, h1);
                        g.drawOval(x2, y2, w, h2);                        
                    }
                } else {
                    if (fdraw) { g.setPaint(fillColor);   g.fill(shape); }
                    if (sdraw) { g.setPaint(strokeColor); g.draw(shape); }
                }
        	}
            
        } else if ( shape instanceof Line2D ) {
            if (sdraw) {
                Line2D l = (Line2D)shape;
                x = (int)(l.getX1()+0.5);
                y = (int)(l.getY1()+0.5);
                w = (int)(l.getX2()+0.5);
                h = (int)(l.getY2()+0.5);
                g.setPaint(strokeColor);
                g.drawLine(x, y, w, h);
            }
        } else {
            if (fdraw) { g.setPaint(fillColor);   g.fill(shape); }
            if (sdraw) { g.setPaint(strokeColor); g.draw(shape); }
        }
        if ( sdraw ) {
            g.setStroke(origStroke);
        }
    }
}
