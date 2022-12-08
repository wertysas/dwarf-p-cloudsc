# -*- coding: utf-8 -*-
from __future__ import annotations
from pydantic import BaseModel, validator
from typing import Any, Dict, Optional, Union, Type


class DataTypes(BaseModel):
    """Specify the datatypes for bool, float and integer fields."""

    bool: Type
    float: Type
    int: Type


class GT4PyConfig(BaseModel):
    """Gather options controlling the compilation and execution of the code generated by GT4Py."""

    backend: str
    backend_opts: Dict[str, Any] = {}
    build_info: Optional[Dict[str, Any]] = None
    device_sync: bool = True
    dtypes: DataTypes = DataTypes(bool=bool, float=float, int=int)
    exec_info: Optional[Dict[str, Any]] = None
    managed: Union[bool, str] = "gt4py"
    rebuild: bool = False
    validate_args: bool = False
    verbose: bool = True

    @validator("exec_info")
    @classmethod
    def set_exec_info(cls, v: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        v = v or {}
        return {**v, "__aggregate_data": True}

    def reset_exec_info(self):
        self.exec_info = {"__aggregate_data": self.exec_info.get("__aggregate_data", True)}

    def with_backend(self, backend: Optional[str]) -> GT4PyConfig:
        args = self.dict()
        if backend is not None:
            args["backend"] = backend
        return GT4PyConfig(**args)

    def with_dtypes(self, dtypes: DataTypes) -> GT4PyConfig:
        args = self.dict()
        args["dtypes"] = dtypes
        return GT4PyConfig(**args)

    def with_validate_args(self, flag: bool) -> GT4PyConfig:
        args = self.dict()
        args["validate_args"] = flag
        return GT4PyConfig(**args)