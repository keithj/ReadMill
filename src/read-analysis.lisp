;;;
;;; Copyright (c) 2010-2011 Genome Research Ltd. All rights reserved.
;;;
;;; This file is part of readmill.
;;;
;;; This program is free software: you can redistribute it and/or modify
;;; it under the terms of the GNU General Public License as published by
;;; the Free Software Foundation, either version 3 of the License, or
;;; (at your option) any later version.
;;;
;;; This program is distributed in the hope that it will b useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;; GNU General Public License for more details.
;;;
;;; You should have received a copy of the GNU General Public License
;;; along with this program.  If not, see <http://www.gnu.org/licenses/>.
;;;

(in-package :uk.ac.sanger.readmill)

(defun write-quality-plot (plot-filespec bam-filespec
                           &key read-group regions index)
  "Writes a mean base quality plot for BAM file BAM-FILESPEC to PLOT-FILESPEC.
Will accept reads of mixed length without warning.

Arguments:

- plot-filespec (pathname designator): The file to which the plot
  will be written. If the file exists, it will be overwritten.
- bam-filespec (pathname designator): The BAM file to be analysed.

Key:

- read-group (string designator): A read group ID to restrict analysis.
- regions (list of region designators): A list of designators for regions
  on reference sequences.
- index (pathname designator): A BAM index, required when regions are
  supplied.

Returns:

- plot-filespec."
  (with-bam (bam (header) (maybe-standard-stream bam-filespec) :regions regions
                 :index index)
    (let ((bam (if read-group
                   (discarding-if (complement (make-rg-p read-group)) bam)
                   bam))
          (title (if read-group
                     (rg-description header read-group)
                     "all read groups")))
      (do* ((aln (next bam) (next bam))
            (rlen 0 (if aln
                        (read-length aln)
                        rlen))
            (sums (make-array rlen :element-type 'fixnum :initial-element 0)
                  (if (> rlen (length sums))
                      (replace (make-array rlen :element-type 'fixnum
                                           :initial-element 0) sums)
                      sums))
            (count 0 (1+ count)))
           ((not (has-more-p bam))
            (plot-per-base plot-filespec (vector-mean sums count)
                           title "Mean base quality"))
        (loop
           with qual = (quality-string aln)
           for i from 0 below rlen
           do (incf (aref sums i) (decode-phred-quality (char qual i)))))))
  plot-filespec)

(defun write-pattern-report (report-filespec pattern-char min-freq
                             bam-filespec &key read-group regions index)
  "Writes a text report describing commonly found patterns of
PATTERN-CHAR in the reads of BAM file denoted by
BAM-FILESPEC. Patterns found MIN-FREQ or more times will be reported.
If READ-GROUP is provided, restricts analysis to those reads only.

Arguments:

- report-filespec (pathname designator): The file to which the report
  will be written. If the file exists, it will be overwritten.
- pattern-char (character): The character whose pattern of
  distribution will be reported.
- min-freq (fixnum): The minimum pattern frequency to be reported.
- bam-filespec (pathname designator): The BAM file to be analysed.

Key:

- read-group (string designator): A read group ID to restrict analysis.
- regions (list of region designators): A list of designators for regions
  on reference sequences.
- index (pathname designator): A BAM index, required when regions are
  supplied.

Returns:

- report-filespec."
  (with-bam (bam (header) (maybe-standard-stream bam-filespec) :regions regions
                 :index index)
    (let ((bam (if read-group
                   (discarding-if (complement (make-rg-p read-group)) bam)
                   bam))
          (title (if read-group
                     (rg-description header read-group)
                     "all read groups")))
      (multiple-value-bind (patterns read-count)
          (base-patterns bam (char-upcase pattern-char) #\.)
        (flet ((report (stream)
                 (format stream "Report for ~a~%" title)
                 (format stream "Total no. reads: ~a~%" read-count)
                 (format stream (txt "~%Most frequent patterns (occurring"
                                  "~a or more times):~%") min-freq)
                 (format stream "~12@a Pattern~%" "Frequency")
                 (dolist (pat (sort patterns #'> :key #'cdr))
                   (when (<= min-freq (cdr pat))
                     (format stream "~12d ~a~%" (cdr pat) (car pat))))))
          (let ((out (maybe-standard-stream report-filespec)))
            (if (streamp out)
                (report out)
                (with-open-file (stream out :direction :output
                                        :if-exists :supersede)
                  (report out)))
            out))))))

(defun base-patterns (bam char mask-char &optional (min-char-freq 0))
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
    (loop
       with pos-table = (make-hash-table :test #'equal)
       while (has-more-p bam)
       for aln = (next bam)
       count aln into read-count
       maximize (read-length aln) into read-length of-type fixnum
       do (let ((seq (seq-string aln)))
            (declare (type simple-base-string seq))
            (loop for i from 0 below (length seq)
               when (char= char (char seq i))
               collect i into pos
               and count i into n
               finally (when (>= n min-char-freq)
                         (if (gethash pos pos-table)
                             (incf (the fixnum (gethash pos pos-table)))
                             (setf (gethash pos pos-table) 1)))))
       finally (return (loop
                          for pos being the hash-keys in pos-table
                          using (hash-value freq)
                          collect (cons (make-mask read-length pos) freq)
                          into summary
                          finally (return (values summary read-count)))))))

(defun rg-description (header read-group)
  (let ((rg (find read-group (header-records (make-sam-header header) :rg)
                  :key (lambda (record)
                         (header-value record :id)) :test #'string=)))
    (check-arguments rg (read-group) "this read-group is not present")
    (format nil "read group ~a ~a~@[ (~a)~]" read-group
            (header-value rg :sm) (header-value rg :pu))))
