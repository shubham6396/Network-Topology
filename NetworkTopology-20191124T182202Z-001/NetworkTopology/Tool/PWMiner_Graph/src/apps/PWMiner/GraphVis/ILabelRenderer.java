package apps.PWMiner.GraphVis;

import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.geom.AffineTransform;
import java.awt.geom.RectangularShape;

import prefuse.Constants;
import prefuse.render.LabelRenderer;
import prefuse.util.ColorLib;
import prefuse.util.FontLib;
import prefuse.util.GraphicsLib;
import prefuse.visual.VisualItem;

public class ILabelRenderer extends LabelRenderer {
	public ILabelRenderer(String textField) {
        super(textField);
    }
    /** Transform used to scale and position images */
    AffineTransform m_transform = new AffineTransform();
    /**
     * @see prefuse.render.Renderer#render(java.awt.Graphics2D, prefuse.visual.VisualItem)
     */
    public void render(Graphics2D g, VisualItem item) {
        RectangularShape shape = (RectangularShape)getShape(item);
        if ( shape == null ) return;
        
        // fill the shape, if requested
        int type = getRenderType(item);
        if ( type==RENDER_TYPE_FILL || type==RENDER_TYPE_DRAW_AND_FILL )
            IGraphicsLib.paint(g, item, shape, getStroke(item), RENDER_TYPE_FILL, item.getInt(CommonDef.TYPE));

        // now render the image and text
        String text = m_text;
        Image  img  = getImage(item);
        
        if ( text == null && img == null )
            return;
                        
        double size = item.getSize();
        boolean useInt = 1.5 > Math.max(g.getTransform().getScaleX(),
                                        g.getTransform().getScaleY());
        double x = shape.getMinX() + size*m_horizBorder;
        double y = shape.getMinY() + size*m_vertBorder;
        
        // render image
        if ( img != null ) {            
            double w = size * img.getWidth(null);
            double h = size * img.getHeight(null);
            double ix=x, iy=y;
            
            // determine one co-ordinate based on the image position
            switch ( m_imagePos ) {
            case Constants.LEFT:
                x += w + size*m_imageMargin;
                break;
            case Constants.RIGHT:
                ix = shape.getMaxX() - size*m_horizBorder - w;
                break;
            case Constants.TOP:
                y += h + size*m_imageMargin;
                break;
            case Constants.BOTTOM:
                iy = shape.getMaxY() - size*m_vertBorder - h;
                break;
            default:
                throw new IllegalStateException(
                        "Unrecognized image alignment setting.");
            }
            
            // determine the other coordinate based on image alignment
            switch ( m_imagePos ) {
            case Constants.LEFT:
            case Constants.RIGHT:
                // need to set image y-coordinate
                switch ( m_vImageAlign ) {
                case Constants.TOP:
                    break;
                case Constants.BOTTOM:
                    iy = shape.getMaxY() - size*m_vertBorder - h;
                    break;
                case Constants.CENTER:
                    iy = shape.getCenterY() - h/2;
                    break;
                }
                break;
            case Constants.TOP:
            case Constants.BOTTOM:
                // need to set image x-coordinate
                switch ( m_hImageAlign ) {
                case Constants.LEFT:
                    break;
                case Constants.RIGHT:
                    ix = shape.getMaxX() - size*m_horizBorder - w;
                    break;
                case Constants.CENTER:
                    ix = shape.getCenterX() - w/2;
                    break;
                }
                break;
            }
            
            if ( useInt && size == 1.0 ) {
                // if possible, use integer precision
                // results in faster, flicker-free image rendering
                g.drawImage(img, (int)ix, (int)iy, null);
            } else {
                m_transform.setTransform(size,0,0,size,ix,iy);
                g.drawImage(img, m_transform, null);
            }
        }
        
        // render text
        int textColor = item.getTextColor();
        if ( text != null && ColorLib.alpha(textColor) > 0 ) {
            g.setPaint(ColorLib.getColor(textColor));
            g.setFont(m_font);
            FontMetrics fm = DEFAULT_GRAPHICS.getFontMetrics(m_font);
            Font i_font = FontLib.getFont(m_font.getName(), Font.ITALIC, size*m_font.getSize());

            FontMetrics fm2 = DEFAULT_GRAPHICS.getFontMetrics(i_font);
            // compute available width
            double tw;
            switch ( m_imagePos ) {
            case Constants.TOP:
            case Constants.BOTTOM:
                tw = shape.getWidth() - 2*size*m_horizBorder;
                break;
            default:
                tw = m_textDim.width;
            }
            
            // compute available height
            double th;
            switch ( m_imagePos ) {
            case Constants.LEFT:
            case Constants.RIGHT:
                th = shape.getHeight() - 2*size*m_vertBorder;
                break;
            default:
                th = m_textDim.height;
            }
            
            // compute starting y-coordinate
            y += fm.getAscent();
            switch ( m_vTextAlign ) {
            case Constants.TOP:
                break;
            case Constants.BOTTOM:
                y += th - m_textDim.height;
                break;
            case Constants.CENTER:
                y += (th - m_textDim.height)/2;
            }
            
            // render each line of text. The even line will be writen in normal way but the odd line will be writen in Italic way
            int lh = fm.getHeight(); // the line height
            int start = 0, end = text.indexOf(m_delim);
            int count =0;
            for ( ; end >= 0; y += lh ) {
                drawString(g, fm, fm2, text.substring(start, end), useInt, x, y, tw, count);
                start = end+1;
                end = text.indexOf(m_delim, start); 
                count++;
            }
            drawString(g, fm, fm2, text.substring(start), useInt, x, y, tw, count);
        }
    
        // draw border
        if (type==RENDER_TYPE_DRAW || type==RENDER_TYPE_DRAW_AND_FILL) {
            GraphicsLib.paint(g,item,shape,getStroke(item),RENDER_TYPE_DRAW);
        }
    }
    
    private final void drawString(Graphics2D g, FontMetrics pfm, FontMetrics ifm, String text,
            boolean useInt, double x, double y, double w, int count)
    {
    	FontMetrics fm = null;
    	if (count % 2 == 1)
    		fm = ifm;
    	else
    		fm = pfm;
    	
        // compute the x-coordinate
        double tx;
        switch ( m_hTextAlign ) {
        case Constants.LEFT:
            tx = x;
            break;
        case Constants.RIGHT:
            tx = x + w - fm.stringWidth(text);
            break;
        case Constants.CENTER:
            tx = x + (w - fm.stringWidth(text)) / 2;
            break;
        default:
            throw new IllegalStateException(
                    "Unrecognized text alignment setting.");
        }
        // use integer precision unless zoomed-in
        // results in more stable drawing
        if ( useInt ) {
            g.drawString(text, (int)tx, (int)y);
        } else {
            g.drawString(text, (float)tx, (float)y);
        }
    }
}
