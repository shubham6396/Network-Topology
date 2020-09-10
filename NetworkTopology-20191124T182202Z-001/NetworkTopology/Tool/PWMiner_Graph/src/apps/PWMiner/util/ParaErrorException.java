package apps.PWMiner.util;

public class ParaErrorException extends Exception {
	public ParaErrorException(){
		super("Parameter Error !");
	}
	
	public ParaErrorException(String msg){
		super(msg);
	}
}
