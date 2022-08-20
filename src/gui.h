
#pragma once

#include <GLFW/glfw3.h>

GLFWwindow* GuiInit();
void GuiTerminate(GLFWwindow* window);
void GuiNewFrame();
void ClearBackGround(GLFWwindow* window);
void Render(GLFWwindow* window);