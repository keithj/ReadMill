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

(defparameter *readmill-commands* (make-hash-table :test #'equal)
  "A hash-table mapping a ReadMill command name to a list containing
a CLI object and a function implementing the desired command.")

(defun get-cli (cmd)
  "Returns the CLI object for command string CMD."
  (first (gethash cmd *readmill-commands*)))

(defun get-fn (cmd)
  "Returns the function for command string CMD."
  (second (gethash cmd *readmill-commands*)))

(defun register-command (cmd cli fn)
  "Registers string CMD as calling FN via CLI. CLI must be a symbol
designating a CLI class."
  (check-arguments (and (stringp cmd) (symbolp cli) (functionp fn))
                   (cmd cli fn)
                   "expected a command string, a cli symbol and a function")
  (setf (gethash cmd *readmill-commands*) (list (make-instance cli) fn)))

(defun print-avail-commands (&optional (stream *error-output*))
  "Prints help for all registered commands to STREAM."
  (let ((commands (sort (loop
                           for cmd being the hash-keys of *readmill-commands*
                           collect cmd) #'string<)))
    (format stream "ReadMill version ~s~%~%" *readmill-version*)
    (write-line "Available commands:" stream)
    (terpri stream)
    (mapc (lambda (cmd)
            (let* ((cli (get-cli cmd))
                   (doc (documentation-of cli)))
              (if doc
                  (help-message cli doc stream)
                  (warn "No help was found for ~a" cmd)))) commands)))

(defun readmill-cli ()
  "Applies the appropriate command line interface."
  (flet ((errmsg (c)
           (princ "Error: " *error-output*)
           (write-line (string-capitalize (format nil "~a" c) :end 1)
                       *error-output*)
           (terpri *error-output*)))
    (with-argv (argv)
      (let* ((cmd (first argv))
             (args (rest argv))
             (cli (get-cli cmd)))
        (handler-case
            (if (or (null cmd) (null cli))
                (error 'unknown-command :cli cli :command cmd)
                (funcall (get-fn cmd)
                         (parse-command-line cli args) argv))
          (unknown-command (condition)
            (errmsg condition)
            (print-avail-commands)
            (quit-lisp :status 2))
          (cli-error (condition)
            (errmsg condition)
            (princ "Usage: " *error-output*)
            (let ((msg (documentation-of cli)))
              (if msg
                  (help-message cli msg *error-output*)
                  (warn "No help was found for ~a~%" cmd)))
            (quit-lisp :status 3)))))))

(define-cli verbosity-mixin ()
  ((verbose "verbose" :required-option nil :value-type t
            :documentation "Print summary of action to standard output.")))

(define-cli input-mixin ()
  ((input "input" :required-option t :value-type 'string
          :documentation "The input file or stream.")))

(define-cli output-mixin ()
  ((output "output" :required-option t :value-type 'string
           :documentation "The output file or stream.")))

(define-cli json-log-mixin ()
  ((json-file "json-file" :required-option nil :value-type 'string
              :documentation "The JSON output file.")))

(define-cli sample-name-mixin ()
  ((sample-name "sample-name" :required-option t :value-type 'string
                :documentation "The sample name.")))

(define-cli read-group-mixin ()
  ((read-group "read-group" :required-option nil :value-type 'string
               :documentation "The restrict results to one read group.")))

(define-cli read-range-mixin ()
  ((read-start "read-start" :required-option nil :value-type 'integer
               :documentation "The zero-based start position in the read.")
   (read-end "read-end" :required-option nil :value-type 'integer
             :documentation #.(txt "The zero-based, half-open end position"
                                   "in the read."))))

(define-cli about-cli (cli)
  ((platform "platform" :required-option nil :value-type t
             :documentation "Reports the Lisp platform details.")
   (version "version" :required-option nil :value-type t
            :documentation "Reports the software version."))
  (:documentation "about [--platform] [--version]"))

(define-cli quality-plot-cli (cli read-group-mixin input-mixin)
  ((plot-file "plot-file" :required-option t :value-type 'string
              :documentation "The plot file."))
  (:documentation "quality-plot --input <name> [--read-group <id>]
--plot-file <filename>"))

(define-cli pattern-report-cli (cli read-group-mixin input-mixin)
  ((report-file "report-file" :required-option t :value-type 'string
                :documentation "The plot file.")
   (pattern-char "pattern-char" :required-option t :value-type 'character
                 :documentation "The pattern character.")
   (min-freq "min-freq" :required-option t :value-type 'integer
             :documentation "The minimum pattern frequency to be reported."))
  (:documentation "pattern-report --input <name> [--read-group <id>]
--report-file <filename> --pattern-char <character> --min-freq <frequency>"))

(define-cli read-filter-cli (cli read-range-mixin
                             input-mixin output-mixin json-log-mixin)
  ((min-quality "min-quality" :required-option nil :value-type 'integer
                :documentation "The minimum quality threshold.")
   (queries "queries" :required-option nil :value-type 'string-list
            :documentation "Subsequences to search for.")
   (orphans "orphans" :required-option nil :value-type t
            :documentation "Include orphans in output."))
  (:documentation "read-filter --input <name> --output <name>
[--read-start <integer>] [--read-end <integer>]
[--min-quality <integer>] [--queries <seq1,seq2 ... seqn>]
[--orphans] [--json-log <filename>]"))

(define-cli split-bam-cli (cli input-mixin output-mixin)
  ((separator "separator" :required-option nil :value-type 'string
              :documentation "The output filename separator string.")
   (max-size "max-size" :required-option t :value-type 'integer
             :documentation "The maximum number of reads per output file."))
  (:documentation "split-bam --input <name> --output <name>
--max-size <integer> [--separator <string>]"))

(register-command "about" 'about-cli #'about)
(register-command "pattern-report" 'pattern-report-cli #'pattern-report)
(register-command "quality-plot" 'quality-plot-cli #'quality-plot)
(register-command "read-filter" 'read-filter-cli #'read-filter)
(register-command "split-bam" 'split-bam-cli #'split-bam)
