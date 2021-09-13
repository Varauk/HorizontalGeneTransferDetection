package objects;

import java.util.ArrayList;

/**
 * Represents a square matrix storing generic elements. Supported matrix formats are (non strict)
 * upper and lower triangular matrices (that is including diagonal elements) in which case only
 * those elements are stored and full matrices containing all n*n elements.
 */
@SuppressWarnings("serial")
public class Matrix<E> extends ArrayList<ArrayList<E>> {

    /**
     * Format of the matrix (i.e. upper triangular, lower triangular or full matrix).
     */
    private MatrixFormat format;

    public Matrix(MatrixFormat format) {
        this.format = format;
    }

    /**
     * Changes the format of the matrix holding the relation.
     *
     * <ul>
     * <li>If current format is BOTH the upper or lower triangle is extracted</li>
     * <li>If current format is UPPER or LOWER and target format is LOWER or UPPER the matrix is
     * transposed</li>
     * <li>If current format is UPPER or LOWER and target format is BOTH the given triangle is
     * reflected and added to the matrix</li>
     * </ul>
     *
     * @param outputFormat New format for the matrix
     */
    public void changeFormat(MatrixFormat outputFormat) {
        // checker whether there is anything to be done
        if (this.format == outputFormat) {
            return;
        }

        switch (this.format) {
            case BOTH:
                // if this matrix is BOTH just delete the correct triangle
                int startOfDelInterval = outputFormat == MatrixFormat.LOWER ? 1 : 0;
                int endOfDelInterval = outputFormat == MatrixFormat.LOWER ? size() : 0;

                for (ArrayList<E> row : this) {
                    row.subList(startOfDelInterval, endOfDelInterval).clear();
                    if (outputFormat == MatrixFormat.LOWER) {
                        startOfDelInterval++;
                    } else {
                        endOfDelInterval++;
                    }
                }
                break;
            case LOWER:
                for (int rowIdx = 0; rowIdx < size(); rowIdx++) {
                    if (outputFormat == MatrixFormat.UPPER) {
                        ArrayList<E> newRow = new ArrayList<E>();
                        for (int colIdx = rowIdx; colIdx < size(); colIdx++) {
                            newRow.add(get(colIdx).get(rowIdx));
                        }
                        set(rowIdx, newRow);
                    } else {   // outputFormat == MatrixFormat.BOTH
                        for (int colIdx = rowIdx + 1; colIdx < size(); colIdx++) {
                            get(rowIdx).add(get(colIdx).get(rowIdx));
                        }
                    }
                }
                break;
            case UPPER:
                for (int rowIdx = size() - 1; rowIdx >= 0; rowIdx--) {
                    if (outputFormat == MatrixFormat.LOWER) {
                        ArrayList<E> newRow = new ArrayList<E>();
                        for (int colIdx = 0; colIdx <= rowIdx; colIdx++) {
                            newRow.add(get(colIdx).get(rowIdx - colIdx));
                        }
                        set(rowIdx, newRow);
                    } else { // outputFormat == MatrixFormat.BOTH
                        for (int colIdx = rowIdx - 1; colIdx >= 0; colIdx--) {
                            get(rowIdx).add(0, get(colIdx).get(rowIdx - colIdx));
                        }
                    }
                }
                break;
            default:
                // should never arise
        }

        // set new format
        this.format = outputFormat;
    }

    public MatrixFormat getFormat() {
        return this.format;
    }

    /**
     * Represents the different forms the matrix holding a gene relation might have.
     */
    public enum MatrixFormat {
        /**
         * Lower triangular matrix.
         */
        LOWER,
        /**
         * Upper triangular Matrix.
         */
        UPPER,
        /**
         * Full matrix, i.e. containing upper and lower triangle.
         */
        BOTH;
    }
}
