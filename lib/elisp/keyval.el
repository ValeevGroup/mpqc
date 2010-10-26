;;
;; keyval.el: Mode for MPQC input files
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
(defvar keyval-mode-hook nil
  "*List of hook functions run by `keyval-mode' (see `run-hooks').")

(defvar keyval-mode-map
  (let ((map (make-sparse-keymap)))
    ;(define-key map " " 'scroll-up)
    ;(define-key map "\^?" 'scroll-down)
    map)
  "Keymap for KeyVal input buffers.")

(defun keyval-mode ()
  "Major mode for KeyVal input files."
  (interactive)
  (fundamental-mode)
  (use-local-map keyval-mode-map)
  (setq major-mode 'keyval-mode
	mode-name "KeyVal")
  (keyval-setup)
  (make-local-variable 'font-lock-defaults)
  (setq font-lock-defaults '(keyval-font-lock-keywords t t))
  (run-hooks 'keyval-mode-hook))

(defun keyval-setup ())

(defvar keyval-key-face 'keyval-key-face
  "Face for keys in KeyVal file.")
(make-face keyval-key-face)
(make-face-bold keyval-key-face)
(set-face-foreground keyval-key-face "Cyan")

(defvar keyval-classtype-face 'keyval-classtype-face
  "Face for names of classes in KeyVal file.")
(make-face keyval-classtype-face)
(set-face-foreground keyval-classtype-face "Green")

(defvar keyval-reference-face 'keyval-reference-face
  "Face for references in KeyVal file.")
(make-face keyval-reference-face)
(set-face-foreground keyval-reference-face "Orange")
(set-face-underline-p keyval-reference-face t)

(defvar keyval-font-lock-keywords
  '(("%.*" . font-lock-comment-face)
    ("<.*>" . keyval-classtype-face)
    ("\"[^\"\n]+\"" . font-lock-string-face)
    ("$[A-Za-z0-9_\.:*+-/]*" . keyval-reference-face)
    ("{ *\\([A-Za-z0-9_\.*+-/ ]*\\>\\) *} *=" (1 keyval-key-face))
    ("\\([A-Za-z0-9_\.*+-/]*\\>\\) *=" (1 keyval-key-face))
    )
  "Default expressions to highlight in KeyVal mode.")
