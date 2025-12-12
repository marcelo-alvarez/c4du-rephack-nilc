"""Preprocessing module for ACT maps: inpainting, filtering, beam corrections, and color corrections."""

from .color_corrections import (
    load_passband,
    compute_color_correction,
    apply_color_correction
)

__all__ = [
    'load_passband',
    'compute_color_correction',
    'apply_color_correction'
]
