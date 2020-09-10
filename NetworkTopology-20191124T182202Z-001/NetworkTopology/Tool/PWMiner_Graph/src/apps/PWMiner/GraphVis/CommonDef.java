package apps.PWMiner.GraphVis;

import javax.swing.ImageIcon;

public class CommonDef {
    
    public static final String LABEL = "label";
    public static final String TYPE = "type";
    
	public static int ISPATTERN = 1;	
	public static int ISTEXT = 2;	
	public static int ISMAPPING = 3;		
	public static int ISNULL = 0;
	public static String PATTERN = "pattern";	
	public static String TEXT = "text";
	public static String MAPPING = "mapping";	
	
	public static String SAVE = "Submit";
	
	public static boolean ISAPPLET = true;
	public static int iAPPLET = 1;
	public static int iAPPLICATION = 2;

	public CommonDef(){
		
	}
	public static ImageIcon createImageIcon(String filename, String description) {
		String path = System.getProperty("user.dir") + "/.." + "/resources/images/" + filename;
		return new ImageIcon(path, description); 
	}
}
