package org.boscaiolo.workout;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import org.junit.Test;
import org.junit.runners.JUnit4;
import org.junit.runner.RunWith;

@RunWith(JUnit4.class)
public class IterTest {

    @Test
    public void test() throws FileNotFoundException {
        BufferedReader reader = new BufferedReader(new FileReader("src/test/resources/sss.txt"));
        IntIter intIter = new IntIter(reader);
        while(intIter.hasNext()){
            System.out.println(intIter.next());
        }
    }
}

