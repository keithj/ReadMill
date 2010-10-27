;;;
;;; Copyright (C) 2010 Genome Research Ltd. All rights reserved.
;;;
;;; This file is part of readmill.
;;;
;;; This program is free software: you can redistribute it and/or modify
;;; it under the terms of the GNU General Public License as published by
;;; the Free Software Foundation, either version 3 of the License, or
;;; (at your option) any later version.
;;;
;;; This program is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;; GNU General Public License for more details.
;;;
;;; You should have received a copy of the GNU General Public License
;;; along with this program.  If not, see <http://www.gnu.org/licenses/>.
;;;

(in-package :uk.ac.sanger.readmill)

(defun write-quality-plot (plot-filespec sample-name bam-filespec)
  "Writes a mean base quality plot for BAM file BAM-FILESPEC to
PLOT-FILESPEC using the plot title SAMPLE-NAME. Will accept reads of
mixed length without warning."
  (with-bgzf (bam bam-filespec)
    (do* ((aln (read-alignment bam) (read-alignment bam))
          (rlen 0 (if aln
                      (read-length aln)
                      rlen))
          (sums (make-array rlen :element-type 'fixnum :initial-element 0)
                (if (> rlen (length sums))
                    (replace (make-array rlen :element-type 'fixnum
                                         :initial-element 0) sums)
                    sums))
          (count 0 (1+ count)))
         ((null aln) (plot-per-base plot-filespec (vector-mean sums count)
                                    sample-name "Mean base quality"))
      (loop
         with qual = (quality-string aln)
         for i from 0 below rlen
         do (incf (aref sums i) (decode-phred-quality (char qual i)))))))

(defun write-pattern-report (report-filespec sample-name pattern-char min-freq
                             bam-filespec)
  "Writes a text report describing commonly found patterns of
PATTERN-CHAR in the reads of BAM file denoted by
BAM-FILESPEC. Patterns found MIN-FREQ or more times will be reported.

Arguments:

- report-filespec (pathname designator): The file to which the report
will be written. If the file exists, it will be overwritten.
- sample-name (string designator): The sample name.
- pattern-char (character): The character whose pattern of
distribution will be reported.
- min-freq (fixnum): The minimum pattern frequency to be reported.
- bam-filespec (pathname designator): The BAM file to be analysed.

Returns:

- T."
  (multiple-value-bind (patterns read-count)
      (base-patterns (char-upcase pattern-char) #\. bam-filespec)
    (with-open-file (out report-filespec :direction :output
                         :if-exists :supersede)
      (format out "Report for sample ~a ~a~%" sample-name
              (parse-namestring bam-filespec))
      (format out "Total no. reads: ~a~%" read-count)
      (format out (txt "~%Most frequent patterns (occurring"
                       "~a or more times):~%") min-freq)
      (format out "~12@a Pattern~%" "Frequency")
      (dolist (pat (sort patterns #'> :key #'cdr))
        (when (<= min-freq (cdr pat))
          (format out "~12d ~a~%" (cdr pat) (car pat))))))
  t)

(defun base-patterns (char mask-char bam-filespec &optional (min-char-freq 1))
  "Returns an alist mapping strings to integers describing the most
common patterns of CHAR in the read sequences from BAM-FILESPEC and
their frequencies. CHAR must occur at least MIN-CHAR-FREQ times in the
sequence to be reported. All other characters in the sequences are
masked with MASK-CHAR. Returns two values; the alist and the total
number of reads examined."
  (declare (optimize (speed 3) (safety 0)))
  (declare (type fixnum min-char-freq))
  (flet ((make-mask (len positions)
           (let ((mask (make-array len :element-type 'base-char
                                  :initial-element mask-char)))
             (loop
                for i in positions
                do (setf (char mask i) char)
                finally (return mask)))))
    (with-bgzf (bam bam-filespec)
      (read-bam-meta bam)
      (loop
         with pos-table = (make-hash-table :test #'equal)
         for aln = (read-alignment bam)
         while aln
         count aln into read-count
         maximize (read-length aln) into read-length of-type fixnum
         do (let ((seq (seq-string aln)))
              (declare (type simple-base-string seq))
              (loop for i from 0 below (length seq)
                 when (char= char (char seq i))
                 count i into n
                 collect i into pos
                 finally (when (>= n min-char-freq)
                           (if (gethash pos pos-table)
                               (incf (the fixnum (gethash pos pos-table)))
                               (setf (gethash pos pos-table) 1)))))
         finally (return
                   (loop
                      for pos being the hash-keys in pos-table
                      using (hash-value freq)
                      collect (cons (make-mask read-length pos) freq)
                      into summary
                      finally (return (values summary read-count))))))))
