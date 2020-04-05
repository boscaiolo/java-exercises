package org.boscaiolo.workout;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class IntIter implements Iterator<Integer> {
                
    private BufferedReader bufferedReader;

    /** The current line. */
    private String cachedLine;

    /** A flag indicating if the iterator has been fully read. */
    private boolean finished = false;

    public IntIter(BufferedReader bufferedReader) {
        this.bufferedReader = bufferedReader;
    }

    @Override
    public boolean hasNext() {
        if (cachedLine != null) {
            return true;
        } else if (finished) {
            return false;
        } else {
            try {
                while (true) {
                    String line = bufferedReader.readLine();
                    if (line == null) {
                        finished = true;
                        return false;
                    } else if (isValidLine(line)) {
                        cachedLine = line;
                        return true;
                    }
                }
            } catch (IOException ioe) {
                return false;
            }
        }
    }

    protected boolean isValidLine(String line) {
        try {
            Integer.parseInt(line);
        } catch (NumberFormatException e) {
            return false;
        }
        return true;
    }

    @Override
    public Integer next() {
        if (!hasNext()) {
            throw new NoSuchElementException("No more lines");
        }
        int current = Integer.parseInt(cachedLine);
        cachedLine = null;
        return current;
    }

}