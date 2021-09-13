package objects;

@SuppressWarnings("serial")
public class NexusFormatException extends Exception {

  public NexusFormatException() {
  }

  public NexusFormatException(String message) {
    super(message);
  }

  public NexusFormatException(String message, Throwable cause) {
    super(message, cause);
  }

  public NexusFormatException(Throwable cause) {
    super(cause);
  }
}
