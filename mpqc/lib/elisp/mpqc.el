;;
;; mpqc.el: mode for editing MPQC output files
;;
;; Copyright (C) 1996 Limit Point Systems, Inc.
;;
;; Author: Curtis Janssen <cljanss@ca.sandia.gov>
;; Maintainer: SNL
;;
;; This file is part of MPQC.
;;
;; MPQC is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation; either version 2, or (at your option)
;; any later version.
;;
;; MPQC is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.
;;
;; You should have received a copy of the GNU General Public License
;; along with the MPQC; see the file COPYING.  If not, write to
;; the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
;;
;; The U.S. Government is granted a limited license as per AL 91-7.
;;

;;;###autoload
(defvar mpqc-mode-hook nil
  "*List of hook functions run by `mpqc-mode' (see `run-hooks').")

(defvar mpqc-mode-map
  (let ((map (make-sparse-keymap)))
    (define-key map " " 'scroll-up)
    (define-key map "\^?" 'scroll-down)
    map)
  "Keymap for mpqc output buffers.")

(defun mpqc-mode ()
  "Major mode for mpqc output files."
  (interactive)
  (fundamental-mode)
  (use-local-map mpqc-mode-map)
  (setq major-mode 'mpqc-mode
	mode-name "Mpqc")
  (mpqc-setup)
  (make-local-variable 'font-lock-defaults)
  (setq font-lock-defaults '(mpqc-font-lock-keywords t))
  (run-hooks 'mpqc-mode-hook))

(defun mpqc-setup ()
  (toggle-read-only 1)
  )

(defvar mpqc-classtype-face 'mpqc-classtype-face
  "Face for names of classes in MPQC output.")
(make-face mpqc-classtype-face)
(set-face-foreground mpqc-classtype-face "Green")

(defvar mpqc-key-face 'mpqc-key-face
  "Face for keys in MPQC output.")
(make-face mpqc-key-face)
(make-face-bold mpqc-key-face)
(set-face-foreground mpqc-key-face "Cyan")

(defvar mpqc-success-face 'mpqc-success-face
  "Face for usable mpqc output.")
(make-face mpqc-success-face)
(set-face-foreground mpqc-success-face "Green")

(defvar mpqc-coor-face 'mpqc-coor-face
  "Face for names of simple internal coordinates in MPQC output.")
(make-face mpqc-coor-face)
(make-face-bold mpqc-coor-face)
(set-face-foreground mpqc-coor-face "Green")

(defvar mpqc-info-face 'mpqc-info-face
  "Face for informational messages in MPQC output.")
(make-face mpqc-info-face)
(set-face-foreground mpqc-info-face "Orange")

(defvar mpqc-warning-face 'mpqc-warning-face
  "Face for warnings in MPQC output.")
(make-face mpqc-warning-face)
(set-face-foreground mpqc-warning-face "Red")

(defvar mpqc-error-face 'mpqc-error-face
  "Face for errors in MPQC output.")
(make-face mpqc-error-face)
(make-face-bold mpqc-error-face)
(set-face-foreground mpqc-error-face "Red")

(defvar mpqc-plain-face 'mpqc-plain-face
  "Face for plain MPQC output.")
(make-face mpqc-plain-face)
(set-face-foreground mpqc-plain-face "White")

(defvar mpqc-font-lock-keywords
  '(
    (".*have been met.*" . mpqc-success-face)
    (".*iter.*$" . mpqc-plain-face)
    ("<.*>" . mpqc-classtype-face)
    ("\"[^\"\n]+\"" . font-lock-string-face)
    ("\\(.*::.*\\) *=\\(.*\\)" (1 mpqc-info-face) (2 mpqc-success-face))
    ("\\(total scf energy\\) = \\(.*\\)"
     (1 mpqc-info-face) (2 mpqc-success-face))
    ("\\(nuclear repulsion energy\\) = \\(.*\\)"
     (1 mpqc-info-face) (2 mpqc-success-face))
    ("\\(taking step of size\\) \\(.*\\)"
     (1 mpqc-info-face) (2 mpqc-success-face))
    ("\\(Value of the .*\\): \\(.*\\)"
     (1 mpqc-info-face) (2 mpqc-success-face))
    ("\\(\\(HOMO\\|LUMO\\) is\\) \\(.*\\) = \\(.*\\)"
     (1 mpqc-info-face) (3 mpqc-key-face) (4 mpqc-success-face))
    ("\\(\\(Max\\|RMS\\) .*\\):" (1 mpqc-info-face))
    ("{ *\\([A-Za-z0-9_\.*+-/ ]*\\>\\) *} *=" (1 mpqc-key-face))
    ("\\([A-Za-z0-9_\.*+-/]*\\>\\) *=" (1 mpqc-key-face))
    (" no$" . mpqc-warning-face)
    (" yes$" . mpqc-success-face)
    (".*has converged.*" . mpqc-success-face)
    (".*has NOT converged.*" . mpqc-error-face)
    ("DEBUG.*" . mpqc-warning-face)
    ("WARNING.*" . mpqc-warning-face)
    ("NOTICE.*" . mpqc-warning-face)
    ("TORS" . mpqc-coor-face)
    ("BEND" . mpqc-coor-face)
    ("STRE" . mpqc-coor-face)
    ("OUTP" . mpqc-coor-face)
    )
  "Default expressions to highlight in MPQC mode.")
