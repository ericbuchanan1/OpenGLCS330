#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#define STB_IMAGE_IMPLEMENTATION
#include "glm/stb_image.h"

using namespace std;

constexpr float M_PI = 3.14159265358979323846f;

// Global variables
const GLuint SCR_WIDTH = 1800;
const GLuint SCR_HEIGHT = 1600;
GLfloat deltaTime = 0.0f; // Time between current frame and last frame
GLfloat lastFrame = 0.0f; // Time of last frame
GLboolean isPerspective = GL_TRUE;

// Camera class definition
class Camera {
public:
	glm::vec3 Position = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 Front = glm::vec3(0.0f, 0.0f, -1.0f);
	glm::vec3 Up;
	glm::vec3 Right;
	glm::vec3 WorldUp = glm::vec3(0.0f, 1.0f, 0.0f);

	GLfloat Yaw = -90.0f;
	GLfloat Pitch = 0.0f;

	GLfloat MovementSpeed = 2.5f;
	GLfloat MouseSensitivity = 0.1f;
	GLfloat Zoom = 45.0f;

	// Constructor with vectors
	Camera(glm::vec3 position, glm::vec3 up) : Position(position), WorldUp(up) {
		updateCameraVectors();
	}

	// Default constructor
	Camera() {
		updateCameraVectors();
	}

	// Returns the view matrix calculated using Euler Angles and the LookAt Matrix
	glm::mat4 GetViewMatrix() {
		return glm::lookAt(Position, Position + Front, Up);
	}

	// Processes input received from any keyboard-like input system
	void ProcessKeyboard(int direction, float deltaTime) {
		float velocity = MovementSpeed * deltaTime;
		if (direction == GLFW_KEY_W)
			Position += Front * velocity;
		if (direction == GLFW_KEY_S)
			Position -= Front * velocity;
		if (direction == GLFW_KEY_A)
			Position -= Right * velocity;
		if (direction == GLFW_KEY_D)
			Position += Right * velocity;
		if (direction == GLFW_KEY_Q)
			Position -= Up * velocity;
		if (direction == GLFW_KEY_E)
			Position += Up * velocity;
	}

	// Processes input received from a mouse input system
	void ProcessMouseMovement(float xoffset, float yoffset, GLboolean constrainPitch = true) {
		xoffset *= MouseSensitivity;
		yoffset *= MouseSensitivity;

		Yaw += xoffset;
		Pitch += yoffset;

		// Make sure that when pitch is out of bounds, screen doesn't get flipped
		if (constrainPitch) {
			if (Pitch > 89.0f)
				Pitch = 89.0f;
			if (Pitch < -89.0f)
				Pitch = -89.0f;
		}

		// Update Front, Right and Up Vectors using the updated Euler angles
		updateCameraVectors();
	}

	// Processes input received from a mouse scroll-wheel event
	void ProcessMouseScroll(float yoffset) {
		MovementSpeed += yoffset;
		if (MovementSpeed < 1.0f)
			MovementSpeed = 1.0f;
		if (MovementSpeed > 10.0f)
			MovementSpeed = 10.0f;
	}

private:
	// Calculates the front vector from the Camera's (updated) Euler Angles
	void updateCameraVectors() {
		// Calculate the new Front vector
		glm::vec3 front;
		front.x = cos(glm::radians(Yaw)) * cos(glm::radians(Pitch));
		front.y = sin(glm::radians(Pitch));
		front.z = sin(glm::radians(Yaw)) * cos(glm::radians(Pitch));
		Front = glm::normalize(front);
		// Also re-calculate the Right and Up vector
		Right = glm::normalize(glm::cross(Front, WorldUp));
		Up = glm::normalize(glm::cross(Right, Front));
	}
};

// Instantiate a Camera object
Camera camera(glm::vec3(0.0f, 0.0f, 3.0f), glm::vec3(0.0f, 1.0f, 0.0f));

// Forward declarations
void processInput(GLFWwindow* window);
void framebuffer_size_callback(GLFWwindow* window, int width, int height);

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	// Ensure we only handle press events
	if (action == GLFW_PRESS || action == GLFW_REPEAT) {
		camera.ProcessKeyboard(key, deltaTime); // deltaTime should be calculated each frame
	}
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos) {
	static float lastX = SCR_WIDTH / 2.0f;
	static float lastY = SCR_HEIGHT / 2.0f;
	static bool firstMouse = true;

	if (firstMouse) {
		lastX = xpos;
		lastY = ypos;
		firstMouse = false;
	}

	double xoffset = xpos - lastX;
	double yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

	lastX = xpos;
	lastY = ypos;

	camera.ProcessMouseMovement(xoffset, yoffset);
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
	camera.ProcessMouseScroll(yoffset);
}

// vertex shader source code
const char* vertexShaderSource = "#version 330 core\n"
"layout (location = 0) in vec3 aPos;\n"
"layout (location = 1) in vec3 aNormal;\n" // Normal vector input
"layout (location = 2) in vec2 aTexCoord;\n"
"out vec3 Normal;\n"  // Pass normal to fragment shader
"out vec3 FragPos;\n" // Pass fragment position to fragment shader
"out vec2 TexCoord;\n"
"uniform mat4 model;\n"  // Model matrix
"uniform mat4 view;\n"   // View matrix
"uniform mat4 projection;\n"  // Projection matrix
"void main()\n"
"{\n"
"   gl_Position = projection * view * model * vec4(aPos, 1.0);\n"
"   FragPos = vec3(model * vec4(aPos, 1.0));\n" // Transform vertex position to world space
"   Normal = mat3(transpose(inverse(model))) * aNormal;\n" // Transform normals to world space
"   TexCoord = aTexCoord;\n"
"}\0";

// fragment shader source code
const char* fragmentShaderSource = "#version 330 core\n"
"out vec4 FragColor;\n"
"in vec3 Normal;\n"
"in vec3 FragPos;\n"
"in vec2 TexCoord;\n"
"uniform vec3 lightPos;\n"
"uniform vec3 viewPos;\n"
"uniform vec3 lightColor;\n"
"uniform sampler2D texture1;\n"
"uniform sampler2D decalTexture;\n"
"uniform vec2 decalUVStart;\n"
"uniform vec2 decalUVEnd;\n"
"uniform bool useDecal;\n"

"struct Material {\n"
"    vec3 ambient;\n"
"    vec3 diffuse;\n"
"    vec3 specular;\n"
"    float shininess;\n"
"};\n"
"uniform Material material;\n"
"void main()\n"
"{\n"

"    vec3 ambient = material.ambient * lightColor;\n" // Ambient

"    vec3 norm = normalize(Normal);\n" // Diffuse
"    vec3 lightDir = normalize(lightPos - FragPos);\n"
"    float diff = max(dot(norm, lightDir), 0.0);\n"
"    vec3 diffuse = lightColor * diff * material.diffuse;\n"

"    vec3 viewDir = normalize(viewPos - FragPos);\n" // Specular
"    vec3 reflectDir = reflect(-lightDir, norm);\n"
"    float spec = pow(max(dot(viewDir, reflectDir), 0.0), material.shininess);\n"
"    vec3 specular = spec * lightColor * material.specular;\n"

"    vec3 lighting = ambient + diffuse + specular;\n"// Combine results

// Base texture color
"    vec4 baseTexColor = texture(texture1, TexCoord);\n"

// Decal texture color
"    vec4 decalTexColor = vec4(0.0);\n"
"    if (useDecal) {\n"
"        decalTexColor = texture(decalTexture, (TexCoord - decalUVStart) / (decalUVEnd - decalUVStart));\n"
"    }\n"

// Decal blending
"    bool isWithinDecalArea = TexCoord.x >= decalUVStart.x && TexCoord.x <= decalUVEnd.x && \n"
"                             TexCoord.y >= decalUVStart.y && TexCoord.y <= decalUVEnd.y;\n"
"    vec4 finalTexColor = isWithinDecalArea && useDecal ? mix(baseTexColor, decalTexColor, decalTexColor.a) : baseTexColor;\n"

"    FragColor = vec4(lighting, 1.0) * finalTexColor; \n"// Final color combination
"}\n";

struct Material {
	glm::vec3 ambient;
	glm::vec3 diffuse;
	glm::vec3 specular;
	float shininess;
};

const Material material1 = {
	glm::vec3(1.0f, 1.0f, 1.0f), // Ambient
	glm::vec3(1.0f, 1.0f, 1.0f), // Diffuse
	glm::vec3(0.5f, 0.5f, 0.5f),  // Specular
	0.0f                         // Shininess
};

const Material material2 = {
	glm::vec3(1.0f, 1.0f, 1.0f), // Ambient
	glm::vec3(1.0f, 1.0f, 1.0f), // Diffuse
	glm::vec3(0.5f, 0.5f, 0.5f),  // Specular
	128.0f                         // Shininess
};

void SetMaterialProperties(const Material& material, GLuint shaderProgram) {
	glUniform3fv(glGetUniformLocation(shaderProgram, "material.ambient"), 1, glm::value_ptr(material.ambient));
	glUniform3fv(glGetUniformLocation(shaderProgram, "material.diffuse"), 1, glm::value_ptr(material.diffuse));
	glUniform3fv(glGetUniformLocation(shaderProgram, "material.specular"), 1, glm::value_ptr(material.specular));
	glUniform1f(glGetUniformLocation(shaderProgram, "material.shininess"), material.shininess);
}

// Vertices / Indices for Plane - Complete
void generatePlane(vector<float>& vertices, vector<unsigned int>& indices, float width, float depth) {
	vertices.clear();
	indices.clear();

	// Each vertex consists of 8 floats: 3 for position, 3 for normal, 2 for texture coordinates
	vertices = {
		// positions // normals // texture coords
		-width / 2.0f, 0.0f,  depth / 2.0f,  0.0f, 1.0f, 0.0f,  0.0f, 1.0f, // Top-left
		 width / 2.0f, 0.0f,  depth / 2.0f,  0.0f, 1.0f, 0.0f,  1.0f, 1.0f, // Top-right
		 width / 2.0f, 0.0f, -depth / 2.0f,  0.0f, 1.0f, 0.0f,  1.0f, 0.0f, // Bottom-right
		-width / 2.0f, 0.0f, -depth / 2.0f,  0.0f, 1.0f, 0.0f,  0.0f, 0.0f, // Bottom-left
	};

	indices = {
		0, 1, 2, // First Triangle
		0, 2, 3  // Second Triangle
	};
}
void setupMesh(unsigned int& VAO, unsigned int& VBO, unsigned int& EBO, const vector<float>& vertices, const vector<unsigned int>& indices) {
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);

	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

	glBindVertexArray(0);
}
// Vertices / Indices for solid cylinder - Complete
void generateSolidCylinder(vector<GLfloat>& vertices, vector<GLuint>& indices, GLfloat radius, GLfloat height, GLint radialSegments, GLint heightSegments) {
	vertices.clear();
	indices.clear();

	// Generate vertices for cylinder sides
	for (int i = 0; i <= heightSegments; ++i) {
		float y = -0.5f * height + (float)i / heightSegments * height;
		float v = (float)i / heightSegments;

		for (int j = 0; j <= radialSegments; ++j) {
			float theta = 2.0f * M_PI * j / radialSegments;
			float x = radius * cos(theta);
			float z = radius * sin(theta);
			float u = (float)j / radialSegments;

			// Vertex position
			vertices.push_back(x);
			vertices.push_back(y);
			vertices.push_back(z);
			// Normal vector (pointing outwards, radial direction)
			vertices.push_back(cos(theta));
			vertices.push_back(0.0f);
			vertices.push_back(sin(theta));
			// Texture coordinates
			vertices.push_back(u);
			vertices.push_back(v);
		}
	}

	// Generate indices for the sides
	for (int i = 0; i < heightSegments; ++i) {
		for (int j = 0; j < radialSegments; ++j) {
			int first = i * (radialSegments + 1) + j;
			int second = first + radialSegments + 1;

			indices.push_back(first);
			indices.push_back(second);
			indices.push_back(first + 1);

			indices.push_back(second);
			indices.push_back(second + 1);
			indices.push_back(first + 1);
		}
	}

	// Generate vertices for top and bottom caps
	int capStartIndex = vertices.size() / 8;
	for (int j = 0; j <= radialSegments; ++j) {
		float theta = 2.0f * M_PI * j / radialSegments;
		float x = radius * cos(theta);
		float z = radius * sin(theta);
		float u = (cos(theta) + 1.0f) / 2.0f;
		float v = (sin(theta) + 1.0f) / 2.0f;

		// Top cap vertex
		vertices.push_back(x);
		vertices.push_back(0.5f * height);
		vertices.push_back(z);
		vertices.push_back(0.0f);
		vertices.push_back(1.0f);
		vertices.push_back(0.0f);
		vertices.push_back(u);
		vertices.push_back(v);

		// Bottom cap vertex
		vertices.push_back(x);
		vertices.push_back(-0.5f * height);
		vertices.push_back(z);
		vertices.push_back(0.0f);
		vertices.push_back(-1.0f);
		vertices.push_back(0.0f);
		vertices.push_back(u);
		vertices.push_back(v);
	}

	// Generate indices for the caps
	for (int j = 0; j < radialSegments; ++j) {
		int topCenter = capStartIndex + j * 2;
		int topNext = topCenter + (j == radialSegments - 1 ? -(2 * radialSegments - 2) : 2);
		int bottomCenter = topCenter + 1;
		int bottomNext = bottomCenter + (j == radialSegments - 1 ? -(2 * radialSegments - 2) : 2);

		// Top cap
		indices.push_back(topCenter);
		indices.push_back(topNext);
		indices.push_back(capStartIndex + 2 * radialSegments);

		// Bottom cap
		indices.push_back(bottomCenter);
		indices.push_back(capStartIndex + 2 * radialSegments + 1);
		indices.push_back(bottomNext);
	}

	// Add center vertices for the caps at the end
	vertices.push_back(0.0f); // x
	vertices.push_back(0.5f * height); // y
	vertices.push_back(0.0f); // z
	vertices.push_back(0.0f); // Normal X
	vertices.push_back(1.0f); // Normal Y
	vertices.push_back(0.0f); // Normal Z
	vertices.push_back(0.5f); // u
	vertices.push_back(0.5f); // v

	vertices.push_back(0.0f); // x
	vertices.push_back(-0.5f * height); // y
	vertices.push_back(0.0f); // z
	vertices.push_back(0.0f); // Normal X
	vertices.push_back(-1.0f);// Normal Y
	vertices.push_back(0.0f); // Normal Z
	vertices.push_back(0.5f); // u
	vertices.push_back(0.5f); // v
}

// Vertices / Indices for sphere - Complete
void generateSphere(vector<float>& vertices, vector<unsigned int>& indices, float radius, int latitudeSegments, int longitudeSegments) {
	vertices.clear();
	indices.clear();

	// Vertices
	for (int lat = 0; lat <= latitudeSegments; ++lat) {
		float theta = lat * M_PI / latitudeSegments;
		float sinTheta = sin(theta);
		float cosTheta = cos(theta);
		float v = (float)lat / latitudeSegments; // V coordinate

		for (int lon = 0; lon <= longitudeSegments; ++lon) {
			float phi = lon * 2 * M_PI / longitudeSegments;
			float sinPhi = sin(phi);
			float cosPhi = cos(phi);
			float u = (float)lon / longitudeSegments; // U coordinate

			float x = cosPhi * sinTheta;
			float y = cosTheta;
			float z = sinPhi * sinTheta;

			// Vertex position
			vertices.push_back(radius * x);
			vertices.push_back(radius * y);
			vertices.push_back(radius * z);
			// Normal vector (same as position vector for a sphere)
			vertices.push_back(x);
			vertices.push_back(y);
			vertices.push_back(z);
			// Texture coordinates
			vertices.push_back(u);
			vertices.push_back(v);
		}
	}

	// Indices
	for (int lat = 0; lat < latitudeSegments; ++lat) {
		for (int lon = 0; lon < longitudeSegments; ++lon) {
			int first = (lat * (longitudeSegments + 1)) + lon;
			int second = first + longitudeSegments + 1;

			indices.push_back(first);
			indices.push_back(second);
			indices.push_back(first + 1);

			indices.push_back(second);
			indices.push_back(second + 1);
			indices.push_back(first + 1);
		}
	}
}

// Vertices / Indices for cylinder with hole / tube - Complete
void generateCylinderWithHole(vector<GLfloat>& vertices, vector<GLuint>& indices, GLfloat outerRadius, GLfloat innerRadius, GLfloat height, GLint radialSegments, GLint heightSegments) {
	vertices.clear();
	indices.clear();

	// Vertices and texture coordinates for the sides of the cylinder
	for (int i = 0; i <= heightSegments; ++i) {
		float y = -0.5f * height + (float)i / heightSegments * height;
		float v = (float)i / heightSegments; // V coordinate

		for (int j = 0; j <= radialSegments; ++j) {
			float theta = 2.0f * M_PI * j / radialSegments;
			float u = (float)j / radialSegments; // U coordinate
			float sinTheta = sin(theta);
			float cosTheta = cos(theta);

			// Outer vertex
			vertices.push_back(outerRadius * cosTheta); // x
			vertices.push_back(y);                      // y
			vertices.push_back(outerRadius * sinTheta); // z
			// Normal for outer surface
			vertices.push_back(cosTheta);
			vertices.push_back(0.0f);
			vertices.push_back(sinTheta);
			// Texture coordinates
			vertices.push_back(u);
			vertices.push_back(v);

			// Inner vertex
			vertices.push_back(innerRadius * cosTheta); // x
			vertices.push_back(y);                      // y
			vertices.push_back(innerRadius * sinTheta); // z
			// Normal for inner surface (inverted)
			vertices.push_back(-cosTheta);
			vertices.push_back(0.0f);
			vertices.push_back(-sinTheta);
			// Texture coordinates
			vertices.push_back(u);
			vertices.push_back(v);
		}
	}

	// Indices for the sides of the cylinder
	for (int i = 0; i < heightSegments; ++i) {
		for (int j = 0; j < radialSegments; ++j) {
			int current = (i * (radialSegments + 1) + j) * 2;
			int next = current + 2;
			int above = current + (radialSegments + 1) * 2;
			int aboveNext = above + 2;

			// Outer surface
			indices.push_back(current);
			indices.push_back(next);
			indices.push_back(aboveNext);

			indices.push_back(current);
			indices.push_back(aboveNext);
			indices.push_back(above);

			// Inner surface
			indices.push_back(current + 1);
			indices.push_back(aboveNext + 1);
			indices.push_back(next + 1);

			indices.push_back(current + 1);
			indices.push_back(above + 1);
			indices.push_back(aboveNext + 1);
		}
	}

	// Adding vertices and indices for the top and bottom caps
	long long baseIndex = vertices.size() / 8; // Adjusted for 8 components per vertex

	// Top and bottom center points for caps
	float topCenterY = 0.5f * height;
	float bottomCenterY = -0.5f * height;

	// Generate vertices for the top and bottom caps
	for (int j = 0; j <= radialSegments; ++j) {
		float theta = 2.0f * M_PI * j / radialSegments;
		float u = (cos(theta) + 1.0f) / 2.0f;
		float v = (sin(theta) + 1.0f) / 2.0f;
		float cosTheta = cos(theta);
		float sinTheta = sin(theta);

		// Top outer vertex
		vertices.push_back(outerRadius * cosTheta);
		vertices.push_back(topCenterY);
		vertices.push_back(outerRadius * sinTheta);
		vertices.push_back(0.0f); // Normal X
		vertices.push_back(1.0f); // Normal Y (upwards)
		vertices.push_back(0.0f); // Normal Z
		vertices.push_back(u);    // Texture coordinate
		vertices.push_back(v);    // Texture coordinate

		// Bottom outer vertex
		vertices.push_back(outerRadius * cosTheta);
		vertices.push_back(bottomCenterY);
		vertices.push_back(outerRadius * sinTheta);
		vertices.push_back(0.0f); // Normal X
		vertices.push_back(-1.0f); // Normal Y (downwards)
		vertices.push_back(0.0f); // Normal Z
		vertices.push_back(u);    // Texture coordinate
		vertices.push_back(v);    // Texture coordinate

		// Top inner vertex
		vertices.push_back(innerRadius * cosTheta);
		vertices.push_back(topCenterY);
		vertices.push_back(innerRadius * sinTheta);
		vertices.push_back(0.0f); // Normal X
		vertices.push_back(1.0f); // Normal Y (upwards)
		vertices.push_back(0.0f); // Normal Z
		vertices.push_back(u);    // Texture coordinate
		vertices.push_back(v);    // Texture coordinate

		// Bottom inner vertex
		vertices.push_back(innerRadius * cosTheta);
		vertices.push_back(bottomCenterY);
		vertices.push_back(innerRadius * sinTheta);
		vertices.push_back(0.0f); // Normal X
		vertices.push_back(-1.0f); // Normal Y (downwards)
		vertices.push_back(0.0f); // Normal Z
		vertices.push_back(u);    // Texture coordinate
		vertices.push_back(v);    // Texture coordinate
	}

	// Generate indices for top and bottom caps
	for (int j = 0; j < radialSegments; ++j) {
		int topOuter1 = baseIndex + j * 4;
		int topOuter2 = baseIndex + ((j + 1) % radialSegments) * 4;
		int bottomOuter1 = topOuter1 + 1;
		int bottomOuter2 = topOuter2 + 1;
		int topInner1 = topOuter1 + 2;
		int topInner2 = topOuter2 + 2;
		int bottomInner1 = topInner1 + 1;
		int bottomInner2 = topInner2 + 1;

		// Top cap
		indices.push_back(topOuter1);
		indices.push_back(topInner2);
		indices.push_back(topInner1);

		indices.push_back(topOuter1);
		indices.push_back(topOuter2);
		indices.push_back(topInner2);

		// Bottom cap
		indices.push_back(bottomOuter1);
		indices.push_back(bottomInner1);
		indices.push_back(bottomInner2);

		indices.push_back(bottomOuter1);
		indices.push_back(bottomInner2);
		indices.push_back(bottomOuter2);
	}
}

// Function to generate cube vertices and indices
void generateCube(vector<GLfloat>& vertices, vector<GLuint>& indices, GLfloat length, GLfloat width, GLfloat height) {
	GLfloat halfLength = length / 2.0f;
	GLfloat halfWidth = width / 2.0f;
	GLfloat halfHeight = height / 2.0f;

	// Cube vertices with positions, normals, and texture coordinates
	GLfloat cubeVertices[] = {
		// positions // normals // texture coordinates Back face
		-halfLength, -halfHeight, -halfWidth,   0.0f, 0.0f, -1.0f,  0.0f, 0.0f,
		halfLength, -halfHeight, -halfWidth,   0.0f, 0.0f, -1.0f,  1.0f, 0.0f,
		halfLength,  halfHeight, -halfWidth,   0.0f, 0.0f, -1.0f,  1.0f, 1.0f,
		-halfLength,  halfHeight, -halfWidth,   0.0f, 0.0f, -1.0f,  0.0f, 1.0f,
		// Front face
		-halfLength, -halfHeight,  halfWidth,   0.0f, 0.0f, 1.0f,   0.0f, 0.0f,
		halfLength, -halfHeight,  halfWidth,   0.0f, 0.0f, 1.0f,   1.0f, 0.0f,
		halfLength,  halfHeight,  halfWidth,   0.0f, 0.0f, 1.0f,   1.0f, 1.0f,
		-halfLength,  halfHeight,  halfWidth,   0.0f, 0.0f, 1.0f,   0.0f, 1.0f,
		// Left face
		-halfLength, -halfHeight, -halfWidth,  -1.0f, 0.0f, 0.0f,   0.0f, 0.0f,
		-halfLength, -halfHeight,  halfWidth,  -1.0f, 0.0f, 0.0f,   1.0f, 0.0f,
		-halfLength,  halfHeight,  halfWidth,  -1.0f, 0.0f, 0.0f,   1.0f, 1.0f,
		-halfLength,  halfHeight, -halfWidth,  -1.0f, 0.0f, 0.0f,   0.0f, 1.0f,
		// Right face
		halfLength, -halfHeight, -halfWidth,   1.0f, 0.0f, 0.0f,   0.0f, 0.0f,
		halfLength, -halfHeight,  halfWidth,   1.0f, 0.0f, 0.0f,   1.0f, 0.0f,
		halfLength,  halfHeight,  halfWidth,   1.0f, 0.0f, 0.0f,   1.0f, 1.0f,
		halfLength,  halfHeight, -halfWidth,   1.0f, 0.0f, 0.0f,   0.0f, 1.0f,
		// Bottom face
		-halfLength, -halfHeight, -halfWidth,   0.0f, -1.0f, 0.0f,  0.0f, 0.0f,
		halfLength, -halfHeight, -halfWidth,   0.0f, -1.0f, 0.0f,  1.0f, 0.0f,
		halfLength, -halfHeight,  halfWidth,   0.0f, -1.0f, 0.0f,  1.0f, 1.0f,
		-halfLength, -halfHeight,  halfWidth,   0.0f, -1.0f, 0.0f,  0.0f, 1.0f,
		// Top face
		-halfLength,  halfHeight, -halfWidth,   0.0f, 1.0f, 0.0f,   0.0f, 0.0f,
		halfLength,  halfHeight, -halfWidth,   0.0f, 1.0f, 0.0f,   1.0f, 0.0f,
		halfLength,  halfHeight,  halfWidth,   0.0f, 1.0f, 0.0f,   1.0f, 1.0f,
		-halfLength,  halfHeight,  halfWidth,   0.0f, 1.0f, 0.0f,   0.0f, 1.0f,
	};

	vertices.assign(cubeVertices, cubeVertices + sizeof(cubeVertices) / sizeof(cubeVertices[0]));

	// Cube indices for 12 triangles
	unsigned int cubeIndices[] = {
		// Back face
		0, 1, 2,
		0, 2, 3,
		// Front face
		4, 5, 6,
		4, 6, 7,
		// Left face
		8, 9, 10,
		8, 10, 11,
		// Right face
		12, 13, 14,
		12, 14, 15,
		// Bottom face
		16, 17, 18,
		16, 18, 19,
		// Top face
		20, 21, 22,
		20, 22, 23
	};

	indices.assign(cubeIndices, cubeIndices + sizeof(cubeIndices) / sizeof(cubeIndices[0]));
}

GLuint loadTexture(const string& path) {
	GLuint textureID;
	glGenTextures(1, &textureID);
	glBindTexture(GL_TEXTURE_2D, textureID);

	// Set texture parameters
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	// Load and generate the texture
	int width;
	int height;
	int nrChannels;
	unsigned char* data = stbi_load(path.c_str(), &width, &height, &nrChannels, 0);
	if (data) {
		GLenum format = (nrChannels == 4) ? GL_RGBA : GL_RGB;
		glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
		glGenerateMipmap(GL_TEXTURE_2D);
	}
	else {
		cerr << "Failed to load texture at path: " << path << endl;
	}
	stbi_image_free(data);

	return textureID;
}

unsigned int compileShader(unsigned int type, const char* source) {
	unsigned int id = glCreateShader(type);
	glShaderSource(id, 1, &source, nullptr);
	glCompileShader(id);

	// Error handling
	int result;
	glGetShaderiv(id, GL_COMPILE_STATUS, &result);
	if (!result) {
		int length;
		glGetShaderiv(id, GL_INFO_LOG_LENGTH, &length);
		auto* message = (char*)alloca(length * sizeof(char));
		glGetShaderInfoLog(id, length, &length, message);
		cerr << "Failed to compile " << (type == GL_VERTEX_SHADER ? "vertex" : "fragment") << " shader" << endl;
		cerr << message << endl;
		glDeleteShader(id);
		return 0;
	}

	return id;
}

unsigned int createShader(const char* vertexShader, const char* fragmentShader) {
	unsigned int program = glCreateProgram();
	unsigned int vs = compileShader(GL_VERTEX_SHADER, vertexShader);
	unsigned int fs = compileShader(GL_FRAGMENT_SHADER, fragmentShader);

	glAttachShader(program, vs);
	glAttachShader(program, fs);
	glLinkProgram(program);
	glValidateProgram(program);

	// Error handling
	int result;
	glGetProgramiv(program, GL_LINK_STATUS, &result);
	if (!result) {
		int length;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &length);
		auto* message = (char*)alloca(length * sizeof(char));
		glGetProgramInfoLog(program, length, &length, message);
		cerr << "Failed to link shader program" << endl;
		cerr << message << endl;
		glDeleteProgram(program);
		return 0;
	}

	glDeleteShader(vs);
	glDeleteShader(fs);

	return program;
}

int main() {
	// GLFW initialization
	if (!glfwInit()) {
		cerr << "Failed to initialize GLFW\n";
		return -1;
	}

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

	GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Eric Buchanan", nullptr, nullptr);
	if (window == nullptr) {
		cerr << "Failed to create GLFW window" << endl;
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window);
	glfwSetKeyCallback(window, key_callback);
	glfwSetCursorPosCallback(window, mouse_callback);
	glfwSetScrollCallback(window, scroll_callback);
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED); // Hide the cursor and capture it.
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	glewExperimental = GL_TRUE;
	if (glewInit() != GLEW_OK) {
		cerr << "Failed to initialize GLEW" << endl;
		return -1;
	}

	glEnable(GL_DEPTH_TEST);
	// Load textures for each object

	GLuint cubeTexture = loadTexture("grey_texture.jpg");
	GLuint cube2Texture = loadTexture("black_texture.jpg");
	GLuint planeTexture = loadTexture("white_texture.jpg");
	GLuint sphereDecalTexture = loadTexture("pingpongball_texture.jpg");
	GLuint sphereTexture = loadTexture("orange_texture.jpg");
	GLuint solidCylinderTexture = loadTexture("blue_texture.jpg");
	GLuint solidCylinderCapTexture = loadTexture("cap_texture.jpg");
	GLuint cylinderWithHoleTexture = loadTexture("braid_texture.jpg");

	unsigned int shaderProgram = createShader(vertexShaderSource, fragmentShaderSource);
	if (shaderProgram == 0) {
		cerr << "Failed to create shader program" << endl;
		return -1;
	}

	// Set up light properties
	glm::vec3 lightPos(1.0f, 1.0f, 2.0f); // Position of the light source
	glm::vec3 lightColor(0.35f, 0.35f, 0.35f); // Color of Light

	// Get uniform locations for the light properties
	int lightPosLoc = glGetUniformLocation(shaderProgram, "lightPos");
	int viewPosLoc = glGetUniformLocation(shaderProgram, "viewPos");
	int lightColorLoc = glGetUniformLocation(shaderProgram, "lightColor");

	// Use the shader program to set uniform variables
	glUseProgram(shaderProgram);
	glUniform3fv(lightPosLoc, 1, glm::value_ptr(lightPos));
	glUniform3fv(lightColorLoc, 1, glm::value_ptr(lightColor));

	// Cube 1 Generation - grey main section Complete
	vector<float> cubeVertices;
	vector<unsigned int> cubeIndices;
	generateCube(cubeVertices, cubeIndices, 2.5f /* length */, 1.0f /* width */, 0.325f /* height */);

	// Create VAO, VBO, EBO and setup mesh for Cube 1
	unsigned int cube1VAO;
	unsigned int cube1VBO;
	unsigned int cube1EBO;
	setupMesh(cube1VAO, cube1VBO, cube1EBO, cubeVertices, cubeIndices);

	// Set up vertex attributes for Cube 1
	glBindVertexArray(cube1VAO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)nullptr);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
	glEnableVertexAttribArray(2);
	glBindVertexArray(0);

	// Cube 2 Generation - black handle rotated closer to pingpongball Complete
	vector<float> cube2Vertices;
	vector<unsigned int> cube2Indices;
	generateCube(cube2Vertices, cube2Indices, 0.1f /* length */, 0.701f /* width */, 0.101f /* height */);

	// Create VAO, VBO, EBO and setup mesh for Cube 2
	unsigned int cube2VAO;
	unsigned int cube2VBO;
	unsigned int cube2EBO;
	setupMesh(cube2VAO, cube2VBO, cube2EBO, cube2Vertices, cube2Indices);

	// Set up vertex attributes for Cube 2
	glBindVertexArray(cube2VAO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)nullptr);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
	glEnableVertexAttribArray(2);
	glBindVertexArray(0);

	// Cube 3 Generation - black handle closer to wristband Complete
	vector<float> cube3Vertices;
	vector<unsigned int> cube3Indices;
	generateCube(cube3Vertices, cube3Indices, 0.1f /* length */, 0.4f /* width */, 0.1f /* height */);

	// Create VAO, VBO, EBO and setup mesh for Cube 3
	unsigned int cube3VAO;
	unsigned int cube3VBO;
	unsigned int cube3EBO;
	setupMesh(cube3VAO, cube3VBO, cube3EBO, cube3Vertices, cube3Indices);

	// Set up vertex attributes for Cube 3
	glBindVertexArray(cube3VAO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)nullptr);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
	glEnableVertexAttribArray(2);
	glBindVertexArray(0);

	// Cube 4 Generation - black handle away from wristband Complete
	vector<float> cube4Vertices;
	vector<unsigned int> cube4Indices;
	generateCube(cube4Vertices, cube4Indices, 0.1f /* length */, 0.4f /* width */, 0.1f /* height */);

	// Create VAO, VBO, EBO and setup mesh for Cube 4
	unsigned int cube4VAO;
	unsigned int cube4VBO;
	unsigned int cube4EBO;
	setupMesh(cube4VAO, cube4VBO, cube4EBO, cube4Vertices, cube4Indices);

	// Set up vertex attributes for Cube 4
	glBindVertexArray(cube4VAO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)nullptr);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
	glEnableVertexAttribArray(2);
	glBindVertexArray(0);

	// Cube 5 Generation - black edge larger side Complete
	vector<float> cube5Vertices;
	vector<unsigned int> cube5Indices;
	generateCube(cube5Vertices, cube5Indices, 0.65f /* length */, 1.001f /* width */, 0.326f /* height */);

	// Create VAO, VBO, EBO and setup mesh for Cube 5
	unsigned int cube5VAO;
	unsigned int cube5VBO;
	unsigned int cube5EBO;
	setupMesh(cube5VAO, cube5VBO, cube5EBO, cube5Vertices, cube5Indices);

	// Set up vertex attributes for Cube 5
	glBindVertexArray(cube5VAO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)nullptr);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
	glEnableVertexAttribArray(2);
	glBindVertexArray(0);

	// Cube 6 Generation - smaller black end Complete
	vector<float> cube6Vertices;
	vector<unsigned int> cube6Indices;
	generateCube(cube6Vertices, cube6Indices, 0.05f /* length */, 1.001f /* width */, 0.326f /* height */);

	// Create VAO, VBO, EBO and setup mesh for Cube 6
	unsigned int cube6VAO;
	unsigned int cube6VBO;
	unsigned int cube6EBO;
	setupMesh(cube6VAO, cube6VBO, cube6EBO, cube6Vertices, cube6Indices);

	// Set up vertex attributes for Cube 6
	glBindVertexArray(cube6VAO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)nullptr);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
	glEnableVertexAttribArray(2);
	glBindVertexArray(0);

	// Plane Generation - Complete
	vector<float> planeVertices;
	vector<unsigned int> planeIndices;
	generatePlane(planeVertices, planeIndices, 10.0f /* width */, 10.0f /* depth */);

	// Create VAO, VBO, EBO and setup mesh for Plane
	unsigned int planeVAO;
	unsigned int planeVBO;
	unsigned int planeEBO;
	setupMesh(planeVAO, planeVBO, planeEBO, planeVertices, planeIndices);

	// Set up vertex attributes for Plane
	glBindVertexArray(planeVAO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)nullptr);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
	glEnableVertexAttribArray(2);
	glBindVertexArray(0);

	// Sphere generation - Complete
	vector<float> sphereVertices;
	vector<unsigned int> sphereIndices;
	generateSphere(sphereVertices, sphereIndices, 1.0f /* radius */, 40 /* latitudeSegments */, 40 /* longitudeSegments */);

	// Create VAO, VBO, EBO
	unsigned int sphereVAO;
	unsigned int sphereVBO;
	unsigned int sphereEBO;
	setupMesh(sphereVAO, sphereVBO, sphereEBO, sphereVertices, sphereIndices);
	// Set vertex attributes specific to the sphere
	glBindVertexArray(sphereVAO);
	// Position attribute
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)nullptr);
	glEnableVertexAttribArray(0);
	// Normal attribute
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	// Texture coordinate attribute
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
	glEnableVertexAttribArray(2);
	glBindVertexArray(0);

	// Cylinder generation - Cap Complete
	vector<float> solidCylinderVertices;
	vector<unsigned int> solidCylinderIndices;
	generateSolidCylinder(solidCylinderVertices, solidCylinderIndices, 1.0f/* radius */, 0.53f/* height */, 36/* radial segments */, 1/* height segments */);

	// Create VAO, VBO, EBO and setup mesh
	unsigned int solidCylinderVAO;
	unsigned int solidCylinderVBO;
	unsigned int solidCylinderEBO;
	setupMesh(solidCylinderVAO, solidCylinderVBO, solidCylinderEBO, solidCylinderVertices, solidCylinderIndices);

	// Set vertex attributes specific to the solid cylinder
	glBindVertexArray(solidCylinderVAO);
	// Position attribute
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)nullptr);
	glEnableVertexAttribArray(0);
	// Normal attribute
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	// Texture coordinate attribute
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
	glEnableVertexAttribArray(2);
	glBindVertexArray(0);

	// Cylinder 2 generation - lipbalm base Complete
	vector<float> solidCylinder2Vertices;
	vector<unsigned int> solidCylinder2Indices;
	generateSolidCylinder(solidCylinder2Vertices, solidCylinder2Indices, 1.00f/* radius */, 1.10f/* height */, 36/* radial segments */, 1/* height segments */);

	// Create VAO, VBO, EBO and setup mesh
	unsigned int solidCylinder2VAO;
	unsigned int solidCylinder2VBO;
	unsigned int solidCylinder2EBO;
	setupMesh(solidCylinder2VAO, solidCylinder2VBO, solidCylinder2EBO, solidCylinder2Vertices, solidCylinder2Indices);

	// Set vertex attributes specific to the solid cylinder
	glBindVertexArray(solidCylinder2VAO);
	// Position attribute
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)nullptr);
	glEnableVertexAttribArray(0);
	// Normal attribute
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	// Texture coordinate attribute
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
	glEnableVertexAttribArray(2);
	glBindVertexArray(0);

	// Cylinder with hole generation - Complete
	vector<float> cylinderVertices;
	vector<unsigned int> cylinderIndices;
	generateCylinderWithHole(cylinderVertices, cylinderIndices, 2.0f /* outerRadius */, 1.8f /* innerRadius */, 0.6f /* height */, 36 /* radialSegments */, 36 /* heightSegments */);

	// Create VAO, VBO, EBO and setup mesh for Cylinder with Hole
	unsigned int cylinderVAO;
	unsigned int cylinderVBO;
	unsigned int cylinderEBO;
	setupMesh(cylinderVAO, cylinderVBO, cylinderEBO, cylinderVertices, cylinderIndices);

	// Set up vertex attributes for Cylinder with Hole
	glBindVertexArray(cylinderVAO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)nullptr);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
	glEnableVertexAttribArray(2);
	glBindVertexArray(0);

	// Render loop
	while (!glfwWindowShouldClose(window)) {
		double currentFrame = glfwGetTime();
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		// Input
		processInput(window);

		// Render
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Set view and projection matrices
		glm::mat4 view = camera.GetViewMatrix();
		glm::mat4 projection;
		if (isPerspective) {
			projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
		}
		else {
			float orthoWidth = 10.0f;
			float aspectRatio = (float)SCR_WIDTH / (float)SCR_HEIGHT;
			projection = glm::ortho(-orthoWidth * aspectRatio, orthoWidth * aspectRatio, -orthoWidth, orthoWidth, -10.0f, 10.0f);
		}
		// Update camera/view position
		glUniform3fv(viewPosLoc, 1, glm::value_ptr(camera.Position));

		glUseProgram(shaderProgram);
		SetMaterialProperties(material2, shaderProgram);
		int modelLoc = glGetUniformLocation(shaderProgram, "model");
		int viewLoc = glGetUniformLocation(shaderProgram, "view");
		int projLoc = glGetUniformLocation(shaderProgram, "projection");

		// Set view and projection matrices here
		glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
		glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(projection));

		// Render Solid Cylinder 1 Cap of lipbalm
		glm::vec3 rotationAxis = glm::vec3(0.0f, 1.0f, 0.0f); // Rotate around y-axis
		glm::mat4 solidCylinderModel = glm::translate(glm::mat4(1.0f), glm::vec3(2.0f, 0.40f, -2.0f));
		solidCylinderModel = glm::rotate(solidCylinderModel, 34.0f, rotationAxis);
		glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(solidCylinderModel));

		glBindVertexArray(solidCylinderVAO);
		// Calculate counts for base and cap
		int totalVertexCount = solidCylinderIndices.size();
		int capVertexCount = 36 * 6;
		int baseVertexCount = totalVertexCount - capVertexCount;
		// Render the base of the cylinder
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, solidCylinderTexture); // Bind base texture
		glUniform1i(glGetUniformLocation(shaderProgram, "texture1"), 0);
		glDrawElements(GL_TRIANGLES, baseVertexCount, GL_UNSIGNED_INT, nullptr); // Draw the base

		// Render the top cap of the cylinder
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, solidCylinderCapTexture);// Bind cap texture
		glUniform1i(glGetUniformLocation(shaderProgram, "texture1"), 0);
		// Calculate the offset where the cap vertices start
		intptr_t capOffset = baseVertexCount * sizeof(unsigned int); // Offset in bytes
		glDrawElements(GL_TRIANGLES, capVertexCount, GL_UNSIGNED_INT, (void*)capOffset);// Draw the top cap
		glBindVertexArray(0);

		// Render Solid Cylinder 2 Base of lipbalm
		glBindTexture(GL_TEXTURE_2D, solidCylinderTexture);
		glBindVertexArray(solidCylinder2VAO);
		glm::mat4 solidCylinder2Model = glm::translate(glm::mat4(1.0f), glm::vec3(2.0f, -0.45f, -2.0f));
		glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(solidCylinder2Model));
		glDrawElements(GL_TRIANGLES, solidCylinderIndices.size(), GL_UNSIGNED_INT, nullptr);
		SetMaterialProperties(material1, shaderProgram);// Change material

		// Render Cube - Light Grey Base
		glBindTexture(GL_TEXTURE_2D, cubeTexture);
		glBindVertexArray(cube1VAO);
		glm::mat4 cube1Model = glm::translate(glm::mat4(1.0f), glm::vec3(-2.1f, -0.73f, 2.0f));
		glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(cube1Model));
		glDrawElements(GL_TRIANGLES, cubeIndices.size(), GL_UNSIGNED_INT, nullptr);

		// Render Cube 2 - Small Black Handle, Long Side
		glBindTexture(GL_TEXTURE_2D, cube2Texture);
		glBindVertexArray(cube2VAO);
		glm::mat4 cube2Model = glm::translate(glm::mat4(1.0f), glm::vec3(-0.70f, -0.73f, 2.0f));
		glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(cube2Model));
		glDrawElements(GL_TRIANGLES, cube2Indices.size(), GL_UNSIGNED_INT, nullptr);

		// Render Cube 3 - small black handle, short side closer farther from wristband
		glBindTexture(GL_TEXTURE_2D, cube2Texture);
		glBindVertexArray(cube3VAO);
		// Rotation transformation
		float angle3 = glm::radians(90.0f); // Rotate
		auto rotationAxis3 = glm::vec3(0.0f, 1.0f, 0.0f); // Rotate around y-axis
		glm::mat4 rotation3 = glm::rotate(glm::mat4(1.0f), angle3, rotationAxis3);
		glm::mat4 translation3 = glm::translate(glm::mat4(1.0f), glm::vec3(-0.85f, -0.73f, 2.3f));
		// Combined model matrix (rotation applied before translation)
		glm::mat4 cube3Model = translation3 * rotation3;
		glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(cube3Model));
		glDrawElements(GL_TRIANGLES, cube3Indices.size(), GL_UNSIGNED_INT, nullptr);

		// Render Cube 4 - small black handle, short side closer to wristband
		glBindTexture(GL_TEXTURE_2D, cube2Texture);
		glBindVertexArray(cube4VAO);
		// Rotation transformation
		float angle4 = glm::radians(90.0f); // Rotate
		auto rotationAxis4 = glm::vec3(0.0f, 1.0f, 0.0f); // Rotate around y-axis
		glm::mat4 rotation4 = glm::rotate(glm::mat4(1.0f), angle4, rotationAxis4);
		glm::mat4 translation4 = glm::translate(glm::mat4(1.0f), glm::vec3(-0.85f, -0.73f, 1.7f));
		// Combined model matrix (rotation applied before translation)
		glm::mat4 cube4Model = translation4 * rotation4;
		glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(cube4Model));
		glDrawElements(GL_TRIANGLES, cube4Indices.size(), GL_UNSIGNED_INT, nullptr);

		// Render Cube 5 - Black Edge, Larger side
		glBindTexture(GL_TEXTURE_2D, cube2Texture);
		glBindVertexArray(cube5VAO);
		glm::mat4 cube5Model = glm::translate(glm::mat4(1.0f), glm::vec3(-3.2f, -0.73f, 2.0f));
		glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(cube5Model));
		glDrawElements(GL_TRIANGLES, cube5Indices.size(), GL_UNSIGNED_INT, nullptr);

		// Render Cube 6 - Black Edge, Smaller side
		glBindTexture(GL_TEXTURE_2D, cube2Texture);
		glBindVertexArray(cube6VAO);
		glm::mat4 cube6Model = glm::translate(glm::mat4(1.0f), glm::vec3(-0.85f, -0.73f, 2.0f));
		glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(cube6Model));
		glDrawElements(GL_TRIANGLES, cube6Indices.size(), GL_UNSIGNED_INT, nullptr);

		// Render Cylinder with Hole
		SetMaterialProperties(material1, shaderProgram); // Change material
		glBindTexture(GL_TEXTURE_2D, cylinderWithHoleTexture);
		glBindVertexArray(cylinderVAO);
		glm::mat4 cylinderModel = glm::translate(glm::mat4(1.0f), glm::vec3(-2.0f, -0.75f, -2.0f));
		glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(cylinderModel));
		glDrawElements(GL_TRIANGLES, cylinderIndices.size(), GL_UNSIGNED_INT, nullptr);

		// Render Plane
		SetMaterialProperties(material2, shaderProgram); // Change material
		glBindVertexArray(planeVAO);

		// Enable the decal usage
		glUniform1i(glGetUniformLocation(shaderProgram, "useDecal"), GL_TRUE);

		// Activate and bind the base texture
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, solidCylinderTexture);
		glUniform1i(glGetUniformLocation(shaderProgram, "texture1"), 0);

		// Activate and bind the decal texture
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, planeTexture); // Replace with your decal texture
		glUniform1i(glGetUniformLocation(shaderProgram, "decalTexture"), 1);

		// Define decal start and end texture coordinates
		glm::vec2 decalStart(0.01f, 0.01f);
		glm::vec2 decalEnd(0.99f, 0.99);

		// Pass decal texture coordinates to the shader
		glUniform2f(glGetUniformLocation(shaderProgram, "decalUVStart"), decalStart.x, decalStart.y);
		glUniform2f(glGetUniformLocation(shaderProgram, "decalUVEnd"), decalEnd.x, decalEnd.y);

		// Set the model matrix and draw the plane
		glm::mat4 planeModel = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, -1.0f, 0.0f));
		glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(planeModel));
		glDrawElements(GL_TRIANGLES, planeIndices.size(), GL_UNSIGNED_INT, nullptr);

		// Disable the decal usage after drawing the plane
		glUniform1i(glGetUniformLocation(shaderProgram, "useDecal"), GL_FALSE);

		// Render Sphere
		glBindVertexArray(sphereVAO);
		// Enable the decal usage
		glUniform1i(glGetUniformLocation(shaderProgram, "useDecal"), GL_TRUE);

		// Activate and bind the base texture
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, sphereTexture);
		glUniform1i(glGetUniformLocation(shaderProgram, "texture1"), 0);

		// Activate and bind the decal texture
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, sphereDecalTexture);
		glUniform1i(glGetUniformLocation(shaderProgram, "decalTexture"), 1);
		// Define decal start and end texture coordinates
		glm::vec2 decal1Start(0.1f, 0.1f);
		glm::vec2 decal1End(0.5f, 0.5f);

		// Pass decal texture coordinates to the shader
		glUniform2f(glGetUniformLocation(shaderProgram, "decalUVStart"), decal1Start.x, decal1Start.y);
		glUniform2f(glGetUniformLocation(shaderProgram, "decalUVEnd"), decal1End.x, decal1End.y);
		// Set the model matrix and draw the sphere
		glm::vec3 rotationAxis5 = glm::vec3(0.0f, 1.0f, 0.0f); // Rotate around y-axis
		glm::mat4 sphereModel = glm::translate(glm::mat4(1.0f), glm::vec3(2.0f, -0.1f, 2.0f));
		sphereModel = glm::rotate(sphereModel, 33.0f, rotationAxis5);
		glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(sphereModel));
		glDrawElements(GL_TRIANGLES, sphereIndices.size(), GL_UNSIGNED_INT, nullptr);
		// Enable the decal usage
		glUniform1i(glGetUniformLocation(shaderProgram, "useDecal"), GL_FALSE);

		// glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	// de-allocate all resources once they've outlived their purpose:
	glDeleteVertexArrays(1, &cube1VAO);
	glDeleteBuffers(1, &cube1VBO);
	glDeleteBuffers(1, &cube1EBO);
	glDeleteVertexArrays(1, &cube2VAO);
	glDeleteBuffers(1, &cube2VBO);
	glDeleteBuffers(1, &cube2EBO);
	glDeleteVertexArrays(1, &cube3VAO);
	glDeleteBuffers(1, &cube3VBO);
	glDeleteBuffers(1, &cube3EBO);
	glDeleteVertexArrays(1, &cube4VAO);
	glDeleteBuffers(1, &cube4VBO);
	glDeleteBuffers(1, &cube4EBO);
	glDeleteVertexArrays(1, &cube5VAO);
	glDeleteBuffers(1, &cube5VBO);
	glDeleteBuffers(1, &cube5EBO);
	glDeleteVertexArrays(1, &planeVAO);
	glDeleteBuffers(1, &planeVBO);
	glDeleteBuffers(1, &planeEBO);
	glDeleteVertexArrays(1, &sphereVAO);
	glDeleteBuffers(1, &sphereVBO);
	glDeleteBuffers(1, &sphereEBO);
	glDeleteVertexArrays(1, &cylinderVAO);
	glDeleteBuffers(1, &cylinderVBO);
	glDeleteBuffers(1, &cylinderEBO);
	glDeleteVertexArrays(1, &solidCylinderVAO);
	glDeleteBuffers(1, &solidCylinderVBO);
	glDeleteBuffers(1, &solidCylinderEBO);
	glDeleteVertexArrays(1, &solidCylinder2VAO);
	glDeleteBuffers(1, &solidCylinder2VBO);
	glDeleteBuffers(1, &solidCylinder2EBO);
	glDeleteProgram(shaderProgram);

	// glfw: terminate, clearing all previously allocated GLFW resources.
	glfwTerminate();
	return 0;
}

void processInput(GLFWwindow* window) {
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);
	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		camera.ProcessKeyboard(GLFW_KEY_W, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		camera.ProcessKeyboard(GLFW_KEY_S, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		camera.ProcessKeyboard(GLFW_KEY_A, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
		camera.ProcessKeyboard(GLFW_KEY_D, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
		camera.ProcessKeyboard(GLFW_KEY_Q, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
		camera.ProcessKeyboard(GLFW_KEY_E, deltaTime);

	static bool p_key_was_pressed = false;
	if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
		if (!p_key_was_pressed) { // Prevent continuous toggles while holding down the key
			isPerspective = !isPerspective;
			p_key_was_pressed = true;
		}
	}
	else {
		p_key_was_pressed = false;
	}
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
	glViewport(0, 0, width, height);
}