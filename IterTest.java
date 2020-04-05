import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;

public class IterTest {

    
    public static void main(String[] args) throws FileNotFoundException {
        BufferedReader reader = new BufferedReader(new FileReader("sss.txt"));
        IntIter intIter = new IntIter(reader);
        while(intIter.hasNext()){
            System.out.println(intIter.next());
        }
    }
}

